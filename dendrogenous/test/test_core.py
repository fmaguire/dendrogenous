#!/usr/bin/env python

import dendrogenous
import dendrogenous.core as core
import unittest
import os
import shutil
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from dendrogenous.test.base import BaseTestCase

class TestCore(BaseTestCase):
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

        cls.settings = {"genome_list" : ["Escherichia_coli_IAI39.fas",
                                         "Escherichia_coli_O157_H7_str._Sakai.fas",
                                         "Nanoarchaeum_equitans_Kin4-M.fas"],
                        "genome_dir" : os.path.join('dendrogenous',
                                                    'test',
                                                    'resources'),
                        "output_dir" : 'core_test_dir',
                        "binary_path": os.path.join('dendrogenous',
                                                    'dependencies'),
                        "db_config" : {"host": "REDACTED",
                                       "db": "new_proteins",
                                       "user": "orchard",
                                       "passwd": "REDACTED"},
                        "minimum_seqs": 3}

        output_dir = cls.settings['output_dir']
        subdirs = ['sequences',
                   'alignments',
                   'masks',
                   'phylogenies', os.path.join('phylogenies', 'coded_named')]
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            for subdir in subdirs:
                os.mkdir(os.path.join(output_dir, subdir))

    def setUp(self):
         self.test = core.dendrogenous(self.test_record,
                                       self.settings)


    def test_init_clean(self):
        """
        Ensure class init works correctly when there is no pre-existing output
        """
        expected_seed = (">YP_025292_1\n"
                         "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF")

        self.assertEqual(self.test.seed, expected_seed)
        self.assertIs(self.test.state, False)


    def test__blast(self):
        """
        Test blast function returns output
        """
        blast_output = self.test._blast("Escherichia_coli_O157_H7_str._Sakai.fas")

        expected_output = self.parse_file(os.path.join(self.test_resources,
                                                       "expected_core_blastp_output.xml"))

        self.assertEqual(blast_output.split(os.linesep), expected_output)


    @unittest.skip('tests local server')
    def test__parse(self):
        """
        Ensure correct parsing of hits locally when connected to server
        """

        blastp_xml = self.parse_file(os.path.join(self.test_resources,
                                                  "expected_core_blastp_output.xml"))

        blastp_output = "\n".join(blastp_xml)

        parsed_blast = self.test._parse_blast(blastp_output)

        expected_parsed_hit_id = '15829270'

        self.assertEqual(parsed_blast[0].id, expected_parsed_hit_id)


    def test_get_seqs(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory
        Ensure state is correctly updated
        """
        #self.test.get_seqs()
        ##for some reason I can't even inspect these
        #self.fail()
        pass


    def test_get_seqs_insufficient(self):
        """
        Test that dendrogenous.get_seqs() correctly renames parsed blast hit
        output file if there are not sufficient numbers of sequences
        Ensure state is correctly updated
        """
        pass
        #self.fail()

    def test_align(self):
        """

        """
        pass
        #self.fail()



    @classmethod
    def tearDownClass(cls):
        """
        Clean up the outputs generated during testing
        """
        if os.path.exists(cls.settings['output_dir']):
            shutil.rmtree(cls.settings['output_dir'])

if __name__ == '__main__':
    unittest.main()
