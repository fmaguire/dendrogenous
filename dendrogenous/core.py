#!/usr/bin/env python

import dendrogenous as dg
from dendrogenous.settings import Settings as SettingsClass
import dendrogenous.utils
import os
import io
import re
import pymysql
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class Dendrogenous():
    """
    Tree generation class with methods for all the steps required to generate
    phylogenies.

    Relies on a pre-existing directory structure


    Generally the idea of running this class to to instantiate it with
    a seq_record e.g.

    current = dendrogenous(seq)
    current.generate_phylogeny()

    Where generate phylogeny will check the current state of the tree
    and

    """

    def __init__(self, seq_record, settings):
        """
        Initialise the class and truncate the accession to 20-chars
        """
        if not isinstance(settings, SettingsClass):
            raise ValueError("Supplied settings is not Settings class")

        if not isinstance(seq_record, SeqRecord):
            raise ValueError("Supplied sequence is not SeqRecord class")

        self.settings = settings
        self.seq_name = self._reformat_accession(seq_record)

        self.seed = ">{0}\n{1}".format(self.seq_name,
                                       seq_record.seq)

        self.seq_hits = os.path.join(self.settings.dir_paths['blast_hits'],
                                    self.seq_name) + '.fas' # file containing blast hits

        self.aligned_seqs = os.path.join(self.settings.dir_paths['alignment'], self.seq_name + ".aligned")
        self.masked_seqs = os.path.join(self.settings.dir_paths['mask'], self.seq_name + ".mask")
        self.phylogeny = os.path.join(self.settings.dir_paths['tree'], self.seq_name + ".tre")

    @staticmethod
    def _reformat_accession(seq_record):
        """
        Reformat accessions to be <=20 chars long and not contain any special chars
        """

        if len(seq_record.id) > 20:
            short_id = seq_record.id[:20]
        else:
            short_id = seq_record.id

        reformatted_id = re.sub('[|,/,\,.]', '_', short_id)
        return reformatted_id



    def __execute_cmd(self, cmd, input_str=None):
        """
        Execute a command using subprocess using a string as stdin
        and returns the stdout as a string if an input_str is provided
        Otherwise just executes cmd normally
        """
        if input_str:
            osnull = open(os.devnull, 'w')
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=osnull)
            (stdout, stderr) = proc.communicate(input_str.encode())
            osnull.close()

        else:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            (stdout, stderr) = proc.communicate()

        return stdout.decode()


    def _blast(self, genome):
        """
        Blast seed sequence against db using BLASTP
        """
        blast_settings = self.settings.blast_settings

        blastp_path = self.settings.binary_paths['blastp']
        blast_cmd = "{0} -db {1} " \
                    " -evalue {2} -max_target_seqs {3}" \
                    " -outfmt 5".format(blastp_path,
                                        genome,
                                        blast_settings['evalue'],
                                        blast_settings['num_seqs'])

        blast_output = dg.utils.execute_cmd(blast_cmd, input_str=self.seed)

        return blast_output


    def _parse_blast(self, blast_output):
        """
        Parses a blast_output files and queries the
        MySQL DB to recover the appropriate sequences
        """

        # no get with default vals as this is a mandatory part of the
        # settings.json - use them to connect to mysql server
        con = pymysql.connect(**self.settings.db_config)
        cur = con.cursor()

        # Parse the xml output returns an iterator of blast record
        # objects over the xml blastoutput file
        blast_hits = NCBIXML.parse(io.StringIO(blast_output))

        hit_id_set = set()
        for blast_record in blast_hits:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit_id_set.add(alignment.hit_def)

        hit_records = []
        for hit_id in hit_id_set:
            cur.execute("SELECT sequence FROM cider "
                        "WHERE protein_ID='{0}'".format(hit_id))
            sequence = cur.fetchone()
            if sequence is None:
                self.settings.logger.error("blast hit protein_ID not in db: {}".format(str(hit_id)))
                continue

            sequence_record = SeqRecord(\
                Seq(sequence[0], IUPAC.protein), id=hit_id, description='')
            hit_records.append(sequence_record)
        con.close()

        return hit_records


    def get_seqs(self):
        """
        Get similar sequences to seed by blasting genomes and parsing the output
        """
        num_hits = 0
        for genome in self.settings.genomes:
            blast_output = self._blast(genome)
            blast_hits = self._parse_blast(blast_output)
            num_hits += len(blast_hits)
            with open(self.seq_hits, 'a') as out_fh:
                SeqIO.write(blast_hits, out_fh, 'fasta')

        if num_hits < self.settings.minimums['min_seqs']:
            os.rename(self.seq_hits,
                      os.path.join(self.settings.dir_paths['blast_fail'], self.seq_name + ".insufficient_hits"))
            self.settings.logger.warning("Too few blastp hits for alignment: {}".format(self.seq_name))

    def align(self):
        """
        Align input seqs using kalign
        """
        if not os.path.exists(self.seq_hits):
            self.get_seqs()


        align_cmd = "kalign -i {0} -o {1}".format(self.seq_hits,
                                                  self.aligned_seqs)
        dg.utils.execute_cmd(align_cmd)

    def mask(self):
        '''
        Mask input alignment using trimal.
        Method reruns trimal with automated setting if mask is too short
        If mask is still too short it is moved to fail directory
        '''

        if not os.path.exists(self.aligned_seqs):
            self.align()

        mask_cmd = "trimal -in {0} -out {1} -nogaps".format(self.aligned_seqs,
                                                            self.masked_seqs)
        dg.utils.execute_cmd(mask_cmd)

        # check if mask is big enough
        mask_length = dg.utils.mask_check(self.masked_seqs)

        if mask_length < settings.cut_off:
            remask_cmd = "trimal -in {0}" \
                         " -out {1} -automated1".format(self.aligned_seqs,
                                                        self.masked_seqs)
            dg.utils.execute_cmd(remask_cmd)

            new_mask_length = dg.utils.mask_check(self.masked_seqs)
            if mask_length < settings.cut_off:
                self.settings.logger.warning("Too few sites hits after mask: {}".format(self.seq_name))
                os.rename(self.masked_seqs,
                          os.path.join(self.settings.dir_paths['mask_fail'],
                                       self.seq_name + ".mask_too_short"))

    def phylogeny(self):
        '''Generate phylogeny from masked seqs using FastTree2'''

        if not os.path.exists(self.masked_seqs):
            self.mask()

        phylogeny_cmd = "FastTree -bionj -slow"\
                        " -quiet -out {1} {0}".format(self.masked_seqs,
                                                      self.phylogeny)
        dg.utils.execute_cmd(phylogeny_cmd)

    #def name(self):
    #    """

    #    """
    #    raise NotImplemented
    #    #if not os.path.exists(PHYLOGENY):
    #    #    self.phylogeny()

