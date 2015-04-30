#!/usr/bin/env python

import numpy as np

import dendrogenous as dg
import ete2
from dendrogenous.settings import Settings as SettingsClass
import dendrogenous.utils
import os
import io
import re
import pymysql
import subprocess

import numpy
from Bio import Phylo
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
        Initialise the class
        """
        if not isinstance(settings, SettingsClass):
            raise ValueError("Supplied settings is not Settings class")

        if not isinstance(seq_record, SeqRecord):
            raise ValueError("Supplied sequence is not SeqRecord class")

        self.settings = settings
        self.seq_name = seq_record.id
        ## seq_record with correct name
        #self.seq_record = seq_record
        #self.seq_record.id = self.seq_name

        self.seed = ">{0}\n{1}\n".format(self.seq_name,
                                       seq_record.seq)

        self.seq_hits = os.path.join(self.settings.dir_paths['blast_hits'],
                                    self.seq_name) + '.fas' # file containing blast hits

        self.aligned_seqs = os.path.join(self.settings.dir_paths['alignment'], self.seq_name + ".afa")
        self.masked_seqs = os.path.join(self.settings.dir_paths['mask'], self.seq_name + ".mask")
        self.phylogeny = os.path.join(self.settings.dir_paths['tree'], self.seq_name + ".tre")
        self.named_phylogeny = os.path.join(self.settings.dir_paths['named'], self.seq_name + ".named_tre")
        self.settings.logger.info("{}: Class Initialised".format(self.seq_name))

    def _dependency(self, dependency_file, dependency_method):
        """
        Check dependency file exists and if it doesn't run the
        method that generates it
        """
        if not os.path.exists(dependency_file):
            dependency_method()

    def _check_output(self, expected_file):
        """
        Check the specified expected file
        has been created and is not empty
        raising a PipeError if it has not
        """
        if not os.path.exists(expected_file):
            raise dg.utils.PipeError("Expected file does not exist: {}".format(expected_file))
        else:
            if not os.path.getsize(expected_file) > 0:
                maformed_file = os.path.join(self.settings.dir_paths['ERROR'], os.path.basename(expected_file))
                os.rename(expected_file, malformed_file)
                raise dg.utils.PipeError("Expected file is empty: {}".format(expected_file))

    def _mask_check(self):
        """
        Returns the length of the mask in the mask file
        Designed for testing if the automated masking needs rerun with different settings
        """
        with open(self.masked_seqs, 'rU') as mask_fh:
            sample_seq = next(SeqIO.parse(mask_fh, "fasta"))
        return len(sample_seq.seq)

    def _get_species_name(self, leaf, db_cursor):
        """
        Query the open database to convert a leaf name (a protein_ID) to
        the appropriate species name.
        """
        if leaf == self.seq_name:
            taxa_name = "SEED SEQUENCE [{}]".format(leaf)
        else:
            query = "SELECT species FROM cider WHERE protein_ID='{}'".format(leaf)
            db_cursor.execute(query)
            returned_taxa_name = db_cursor.fetchone()
            if returned_taxa_name is None:
                taxa_name = 'UNKNOWN TAXA [{}]'.format(leaf)
                self.settings.logger.warning(\
                        "{0}: NameError | Protein ID ({1}) is missing species information".format(\
                            self.seq_name, leaf))
            else:
                taxa_name = "{0} [{1}]".format(returned_taxa_name[0],
                                               leaf)
        return taxa_name

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

        blast_output = dg.utils.execute_cmd(blast_cmd, input_str=self.seed, output_str=True)

        return blast_output

    def _parse_blast(self, blast_output):
        """
        Parses a blast_output files and queries the
        MySQL DB to recover the appropriate sequences
        """

        # no get with default vals as this is a mandatory part of the
        # settings.json - use them to connect to mysql server
        con = pymysql.connect(**self.settings.dbconfig)
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
                self.settings.logger.warning("Blast hit protein_ID not in db: {}".format(hit_id))
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

        with open(self.seq_hits, 'a') as out_fh:
            for genome in self.settings.genomes:
                blast_output = self._blast(genome)
                blast_hits = self._parse_blast(blast_output)
                num_hits += len(blast_hits)
                SeqIO.write(blast_hits, out_fh, 'fasta')

            # add the actual seq seed to blast hits so it
            # will be in phylogenies
            out_fh.write(self.seed)

        self._check_output(self.seq_hits)

        if num_hits < self.settings.minimums['min_seqs']:
            seq_fail_file = os.path.join(self.settings.dir_paths['blast_fail'], self.seq_name + ".insufficient_hits")
            os.rename(self.seq_hits, seq_fail_file)
            self._check_output(seq_fail_file)
            raise dg.utils.GetSeqFail()

        self.settings.logger.info("{}: BLAST Seqs Created".format(self.seq_name))


    def align(self):
        """
        Align input seqs using kalign
        """
        self._dependency(self.seq_hits, self.get_seqs)

        kalign_path = self.settings.binary_paths['kalign']
        align_cmd = "{0} -i {1} -o {2}".format(kalign_path,
                                               self.seq_hits,
                                               self.aligned_seqs)
        dg.utils.execute_cmd(align_cmd)
        self._check_output(self.aligned_seqs)
        self.settings.logger.info("{}: Alignment Created".format(self.seq_name))

    def mask(self):
        '''
        Mask input alignment using trimal.
        Method reruns trimal with automated setting if mask is too short
        If mask is still too short it is moved to fail directory
        '''
        self._dependency(self.aligned_seqs, self.align)

        trimal_path = self.settings.binary_paths['trimal']
        mask_cmd = "{0} -in {1} -out {2} -nogaps".format(trimal_path,
                                                         self.aligned_seqs,
                                                         self.masked_seqs)
        dg.utils.execute_cmd(mask_cmd, debug=True)

        # check if mask is big enough
        mask_length = self._mask_check()

        # if too few sites rerun with automated setting
        if mask_length < self.settings.minimums['min_sites']:
            remask_cmd = "{0} -in {1}" \
                         " -out {2} -automated1".format(trimal_path,
                                                        self.aligned_seqs,
                                                        self.masked_seqs)
            dg.utils.execute_cmd(remask_cmd)
            mask_length = self._mask_check()

            # if still too short move to fail pile
            if mask_length < self.settings.minimums['min_sites']:
                mask_fail_file = os.path.join(self.settings.dir_paths['mask_fail'],
                                              self.seq_name + ".mask_too_short")

                os.rename(self.masked_seqs, mask_fail_file)
                # status is False if the mask fails
                self._check_output(mask_fail_file)
                raise dg.utils.MaskFail()
            else:
                self._check_output(self.masked_seqs)
        else:
            self._check_output(self.masked_seqs)
        self.settings.logger.info("{}: Mask Created".format(self.seq_name))

    def estimate_phylogeny(self):
        '''Generate phylogeny from masked seqs using FastTree2'''
        self._dependency(self.masked_seqs, self.mask)

        fasttree_path = self.settings.binary_paths['FastTree']
        phylogeny_cmd = "{0} -bionj -slow"\
                        " -quiet -out {1} {2}".format(fasttree_path,
                                                      self.phylogeny,
                                                      self.masked_seqs)
        dg.utils.execute_cmd(phylogeny_cmd)
        self._check_output(self.phylogeny)
        self.settings.logger.info("{}: Phylogeny Created".format(self.seq_name))

    def name_phylogeny(self):
        """
        Name sequences on a phylogeny
        """
        self._dependency(self.phylogeny, self.estimate_phylogeny)

        # parse tree using biopython parser
        parsed_tree = ete2.Tree(self.phylogeny)#, 'newick')

        # generate tree_name dict using database
        con = pymysql.connect(**self.settings.dbconfig)
        cur = con.cursor()
        tree_rename_dict = {re.escape(leaf.name): self._get_species_name(leaf.name, cur) \
                                for leaf in parsed_tree.get_leaves()}
        con.close()

        # read tree in as text for easy search replace
        # this could maybe be done more elegantly using the already
        # ETE parsed tree object
        with open(self.phylogeny, 'r') as phy_fh:
            tree_text = phy_fh.read()

        patterns = re.compile("|".join(tree_rename_dict.keys()))

        renaming = lambda m: tree_rename_dict[re.escape(m.group(0))]
        renamed_tree = patterns.sub(renaming, tree_text)

        with open(self.named_phylogeny, 'w') as named_phy_fh:
            named_phy_fh.write(renamed_tree)

        self._check_output(self.named_phylogeny)
        self.settings.logger.info("{}: Phylogeny Named".format(self.seq_name))


    def build_named_phylogeny(self):
        """
        Runner method for class
        """
        try:
            self.name_phylogeny()
        except dg.utils.GetSeqFail:
            self.settings.logger.warning("{}: SeqFail | too few blastp hits for alignment".format(self.seq_name))
        except dg.utils.MaskFail:
            self.settings.logger.warning("{}: MaskFail | too few sites hits after mask".format(self.seq_name))
        except dg.utils.PipeError as E:
            self.settings.logger.error("!!Error in phylogeny generation for {0}: {1}".format(self.seq_name, E.msg))
        finally:
            return



class TreeParser():
    """
    Parse a folder of unlabelled, or a directory
    hierarchy of labelled, newick phylogenies
    using a supplied set of labels
    """

    def __init__(self, categories):
        """
        categories = {bacteria: (bacterial, taxa, names),
                      host: (alveolate taxa labels),
                      ...}
        """

        self.ncbi_taxa = ete2.ncbi_taxonomy.NCBITaxa()


    def _get_taxonomy(self, node_labels):
        """
        Return a taxonomic lineage from a tree label name
        """

        correct_name = self.ncbi_taxa.get_fuzzy_name_translation(node_label)
        taxid = self.ncbi_taxa.get_name_translator

    def build_training(self):
        """
        Return a matrix and labels in a vector
        """
        matrix = self.build_matrix()
        labels = self.encode_labels()
        return matrix, labels

    def build_matrix(self):
        """
        Generate a matrix from list of phylogeny filepaths
        """

        vectors = []
        for phylogeny in phylogenies:



            vectors.append()













