#!/usr/bin/env python

import dendrogenous
import dendrogenous.core as core
import unittest
import os
import shutil
import sys

print(sys.path)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class TestCore(unittest.TestCase):
    """
    Unit tests for the core dendrogeous class and its methods
    """
    @classmethod
    def setUpClass(cls):
        """
        Generate structures required for testing just dendrogenous core class
        """
        cls.test_record = SeqRecord(\
                   Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF",
                   IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein, small")

        cls.settings = {"genome_list" : ["Escherichia_coli_IAI39",
                                         "Escherichia_coli_O157_H7_str._Sakai",
                                         "Nanoarchaeum_equitans_Kin4-M"],
                        "genome_dir" : os.path.join('dendrogenous',
                                                    'test',
                                                    'resources'),
                        "output_dir" : 'core_test_dir'}

        output_dir = cls.settings['output_dir']
        subdirs = ['sequences', os.path.join('sequences', 'unparsed_blast_output'),
                   'alignments',
                   'masks',
                   'phylogenies', os.path.join('phylogenies', 'coded_named')]
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            for subdir in subdirs:
                os.mkdir(os.path.join(output_dir, subdir))

    def test_init_clean(self):
        """
        Ensure class init works correctly when there is no pre-existing output
        """
        test = core.dendrogenous(self.test_record,
                                 self.settings)

        expected_seed = (">YP_025292_1\n"
                         "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF")

        self.assertEqual(test.seed, expected_seed)
        self.assertIs(test.state, False)

        # TODO: Make sure tests for utils functions init calls are tested


    def test__blast(self):
        """
        Blast
        """
        self.fail()
        # only output correct number of temp outputs to output folder
        # make sure these outputs are correct and fixed
        # mock db?


    def test__parse(self):
        """
        Ensure correct parsing of hits
        """
        self.fail()
        # reads correct number of outputs
        # correctly deletes them
        # has correct sequences in parsed output in right folder


    def test_get_seqs(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory
        Ensure state is correctly updated
        """
        self.fail()

    def test_get_seqs_insufficient(self):
        """
        Test that dendrogenous.get_seqs() correctly renames parsed blast hit
        output file if there are not sufficient numbers of sequences
        Ensure state is correctly updated
        """
        self.fail()

    def test_align(self):
        """

        """

    @classmethod
    def tearDownClass(cls):
        """
        Clean up the outputs generated during testing
        """
        if os.path.exists(cls.settings['output_dir']):
            shutil.rmtree(cls.settings['output_dir'])

if __name__ == '__main__':
    unittest.main()
