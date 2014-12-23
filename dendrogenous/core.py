#!/usr/bin/env python

import dendrogenous.utils as utils
import os
import io
import pymysql
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class dendrogenous():
    """
    Tree generation class with methods for all the steps required to generate
    phylogenies
    """

    def __init__(self, seq_record, settings):
        """
        Initialise the class and truncate the accession to 20-chars
        """

        self.settings = settings

        self.seq_name = utils.reformat_accession(seq_record)
        self.seed = ">{0}\n{1}".format(self.seq_name,
                                       seq_record.seq)


        self.output_dir = settings['output_dir']
        # check precomputed state of class by looking for files matching
        # this name in output folders (globs)
        self.state = self._get_state(self.seq_name, self.output_dir)

    def _execute_cmd(self, cmd, input_str=None):
        """
        Execute a command using subprocess using a string as stdin
        and returns the stdout as a string if an input_str is provided
        Otherwise just executes cmd normally
        """
        if input_str is None:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            (stdout, stderr) = proc.communicate()


        else:
            osnull = open(os.devnull, 'w')
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                    stderr=osnull)
            (stdout, stderr) = proc.communicate(input_str.encode())

            osnull.close()

        return stdout.decode()

    def _reformat_accession(self, seq_record):
        """
        Reformat accessions to be <=20 chars long and not contain any special chars
        """

        if len(seq_record.id) > 20:
            short_id = seq_record.id[:19]
        else:
            short_id = seq_record.id

        reformatted_id = re.sub('[|,/,\,.]', '_', short_id)
        return reformatted_id


    def _get_state(self, seq_name, output_dir):
        """
        Check the completed run state of a given seed name in the output dir
        e.g. is there already a generated phylogeny and so on
        """
        return False


    def _blast(self, genome_name):
        """
        Blast seed sequence against db using BLASTP
        """
        blast_settings = self.settings.get('blast_settings', {})
        num_seqs = blast_settings.get('num_seqs', 1)
        e_value = blast_settings.get('e_value', 1e-5)

        genome = os.path.join(self.settings['genome_dir'], genome_name)
        blastp_path = os.path.join(self.settings['binary_path'], 'blastp')
        blast_cmd = "{0} -db {1} " \
                    " -evalue {2} -max_target_seqs {3}" \
                    " -outfmt 5".format(blastp_path,
                                        genome,
                                        e_value,
                                        num_seqs)

        blast_output = self._execute_cmd(blast_cmd, input_str=self.seed)

        return blast_output


    def _parse_blast(self, blast_output):
        """
        Parses a blast_output files and queries the
        MySQL DB to recover the appropriate sequences
        """

        # no get with default vals as this is a mandatory part of the
        # settings.json - use them to connect to mysql server
        db_config = self.settings['db_config']
        con = pymysql.connect(host=db_config['host'],
                              user=db_config['user'],
                              passwd=db_config['passwd'],
                              db=db_config['db'])
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
                with open('protein_IDs_that_dont_exist', 'a+') as fh:
                   fh.write(str(hit_id) + '\n')
                continue
                #raise Exception("Failed to find"
                #                " protein_ID in DB: {0}".format(hit_id))

            sequence_record = SeqRecord(\
                Seq(sequence[0], IUPAC.protein), id=hit_id, description='')
            hit_records.append(sequence_record)
        con.close()

        return hit_records


    def get_seqs(self):
        """
        Get similar sequences to seed by blasting genomes and parsing the output
        """
        self.output_seqs = os.path.join(self.output_dir, "sequences",
                                        self.seq_name) # file containing blast hits
        num_hits = 0
        for genome in self.settings['genome_list']:
            blast_output = self._blast(genome)
            blast_hits = self._parse_blast(blast_output)
            num_hits += len(blast_hits)
            SeqIO.write(blast_hits, self.output_seqs, 'fasta')

        if num_hits < self.settings['minimum_seqs']:
            os.rename(self.output_seqs,
                      self.output_seqs+ '_too_few_seqs_for_alignment')

        assert False

    def align(self):
        """
        Align input seqs using kalign
        """

        if not os.path.exists(ALIGNMENT_FILE):
            self.seqs()

        align_cmd = "kalign -i {0} -o {1}".format(input_seqs,
                                              output_aligned_seqs)
        subprocess.call(align_cmd)

    def mask(self):
        '''Mask input alignment using trimal'''

        if not os.path.exists(alignment_file):
            self.align()

        mask_cmd = "trimal -in {0} -out {1} -nogaps".format(input_aligned_seqs,
                                                            output_masked_seqs)
        accessories.system_call(mask_cmd)

        mask_length = accessories.mask_check(output_masked_seqs)

        if mask_length < cut_off:
            remask_cmd = "trimal -in {0}" \
                         " -out {1} -automated1".format(input_aligned_seqs,
                                                        output_masked_seqs)
            accessories.system_call(remask_cmd)
            new_mask_length = accessories.mask_check(output_masked_seqs)
            if mask_length < cut_off:
                mask_error = "{0} is" \
                             " too short" \
                             " after masking".format(output_masked_seqs)
                os.rename(output_masked_seqs,
                          output_masked_seqs + "_mask_too_short")

                # write to sqlite db


    def phylogeny(self):
        '''Generate phylogeny from masked seqs using FastTree2'''

        if not os.path.exists(MASK_FILE):
            self.mask()

        phylogeny_cmd = "FastTree -bionj -slow" \
                        " -quiet -out {1} {0}".format(input_masked_seqs,
                                                      output_phylogeny)
        subprocess.call(phylogeny_cmd)


    def rename(self):
        """

        """
        if not os.path.exists(PHYLOGENY):
            self.phylogeny()


