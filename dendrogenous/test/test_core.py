#!/usr/bin/env python

import dendrogenous as dg
import dendrogenous.core


import logging
import unittest
import unittest.mock as mock
import os
import shutil
import sys
import pytest
import pickle

from socket import gethostname

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from dendrogenous.test.base import BaseTestCase

class TestCoreInit(BaseTestCase):
    """
    Unit tests for the core dendrogeous init and reform methods
    """
    def setUp(self):
        #self.test_class = dg.core.Dendrogenous(self.test_record,
        #                                       self.settings)
        self.test_record = SeqRecord(\
                   Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF",
                   IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein, small")


    def init_core():
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        test_class = dg.core.Dendrogenous(test_record,
                                          mock_settings)
        return test_class

    def test_init_clean(self):
        """
        Ensure class init works when provided a seqrecord and settings file
        """
        settings = mock.Mock(dg.settings.Settings)
        settings.output_dir = "testdir"

        test_class = dg.core.Dendrogenous(self.test_record,
                                          settings)

        expected_name = "YP_025292_1"
        expected_seed = (">YP_025292_1\n"
                         "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF")

        self.assertEqual(test_class.seed, expected_seed)
        self.assertEqual(test_class.seq_name, expected_name)
        self.assertEqual(test_class.output_dir, settings.output_dir)


    def test_init_bad_seqrecord(self):
        """
        Ensure ValueError is raised if class is instantiated without valid seqrecord
        """
        settings = mock.Mock(dg.settings.Settings)
        settings.output_dir = "testdir"

        invalid_seq = ""

        with self.assertRaises(ValueError):
            dg.core.Dendrogenous(invalid_seq,
                                 settings)

    def test_init_bad_settings(self):
        """
        Ensure ValueError is raised if class is instantiated without valid seqrecord
        """
        invalid_settings = ""

        with self.assertRaises(ValueError):
            dg.core.Dendrogenous(self.test_record,
                                 invalid_settings)

    def test_reformat_accession_method_for_too_long_accessions(self):
        """
        Test reformat accession works as expected
        """
        too_long = SeqRecord(\
                   Seq("X",
                   IUPAC.protein),
                   id="012345678901234567890123456789",
                   name="foo",
                   description="bar, baz")

        truncated = dg.core.Dendrogenous._reformat_accession(too_long)

        self.assertEqual(len(truncated), 20)
        self.assertEqual(truncated, "01234567890123456789")

    def test_reformat_accession_method_for_problematic_characters(self):
        """
        Test reformat accession works as expected
        """
        bad_char = SeqRecord(\
                   Seq("X",
                   IUPAC.protein),
                   id="|blah|t",
                   name="foo",
                   description="bar, baz")

        fixed_chars = dg.core.Dendrogenous._reformat_accession(bad_char)
        self.assertEqual(fixed_chars, "_blah_t")


class TestCoreGetSeqs(BaseTestCase):
    """
    Test the components of the get_seqs method and the _blast and _parse
    functions it relies on
    """

    def setUp(self):
        #self.test_class = dg.core.Dendrogenous(self.test_record,
        #                                       self.settings)
        self.test_record = SeqRecord(\
                   Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF",
                   IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein, small")


    def test__blast_runs(self):
        """
        Make sure the __blast method correctly runs and returns decoded xml output
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.blast_settings = {'num_seqs': 1,
                                        'evalue': 1e-5}
        genome = os.path.join(self.test_resources, "Escherichia_coli_O157_H7_str._Sakai.fas")

        binary_paths = os.path.join("dendrogenous", "dependencies")
        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        expected_output = self.parse_file(os.path.join(self.test_resources,
                                                       "expected_core_blastp_output.xml"))

        blast_output = test_class._blast(genome)
        self.assertEqual(blast_output.split(os.linesep), expected_output)

    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test__parse_blast(self):
        """
        Ensure correct parsing of hits locally when connected to server
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.db_config = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        blastp_xml = self.parse_file(os.path.join(self.test_resources,
                                                  "expected_core_blastp_output.xml"))

        blastp_output = "\n".join(blastp_xml)

        parsed_blast = test_class._parse_blast(blastp_output)

        expected_parsed_hit_id = '15829270'
        expected_parsed_seq = Seq('MLNTC', IUPAC.IUPACProtein())

        self.assertEqual(parsed_blast[0].id, expected_parsed_hit_id)
        self.assertEqual(parsed_blast[0].seq[:5], expected_parsed_seq)

    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test__parse_blast_broken(self):
        """
        Ensure correct parsing of hits locally when connected to server
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        mock_settings.logger = logging.getLogger("test")

        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.db_config = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        blastp_xml = self.parse_file(os.path.join(self.test_resources,
                                                  "broken_blastp_output.xml"))
        blastp_output = "\n".join(blastp_xml)

        parsed_blast = test_class._parse_blast(blastp_output)
        expected_output = []

        self.assertEqual(parsed_blast, expected_output)
        # test logging?
        #self.assertIs(test_class.logger.error, "foobar")


    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_get_seqs(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory - integration test with _blast and _parse
        """

        #configure all the dependencies in a mock settings object
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        mock_settings.dir_paths = {'blast_hits': os.path.join("testdir", 'blast_out')}
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.minimums = {'min_seqs': 3}
        mock_settings.blast_settings = {'num_seqs': 2,
                                        'evalue': 5}
        mock_settings.genomes= [os.path.join(self.test_resources, "Escherichia_coli_O157_H7_str._Sakai.fas"),
                                os.path.join(self.test_resources, "Escherichia_coli_IAI39.fas"),
                                os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")]
        mock_settings.logger = logging.getLogger("test")

        # make output dir that is normally done by runner
        os.mkdir('testdir')
        os.mkdir(mock_settings.dir_paths['blast_hits'])

        # load db settings from secret pickle file
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.db_config = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        test_class.get_seqs()

        expected_output_file = os.path.join('testdir', 'blast_out', 'YP_025292_1.fas')
        self.assertTrue(os.path.exists(expected_output_file))

        with open(expected_output_file, 'r') as out_fh:
            seqs = list(SeqIO.parse(out_fh, 'fasta'))
            print(seqs)
            self.assertEqual(len(seqs), 5)


    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_get_seqs_fails_correctly(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory - integration test with _blast and _parse
        Test it correctly identifies too few hits and moves file to insufficient hits dir
        """

        #configure all the dependencies in a mock settings object
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.output_dir = "testdir"
        mock_settings.dir_paths = {'blast_hits': os.path.join("testdir", 'blast_out'),
                                   'blast_fail': os.path.join("testdir", "blast_out", "insufficient")}
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.minimums = {'min_seqs': 10}
        mock_settings.blast_settings = {'num_seqs': 2,
                                        'evalue': 5}
        mock_settings.genomes= [os.path.join(self.test_resources, "Escherichia_coli_O157_H7_str._Sakai.fas"),
                                os.path.join(self.test_resources, "Escherichia_coli_IAI39.fas"),
                                os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")]
        mock_settings.logger = logging.getLogger("test")

        # make output dir that is normally done by runner
        os.mkdir('testdir')
        os.mkdir(mock_settings.dir_paths['blast_hits'])
        os.mkdir(mock_settings.dir_paths['blast_fail'])

        # load db settings from secret pickle file
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.db_config = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)


        test_class.get_seqs()

        expected_output_file = os.path.join('testdir', 'blast_out', 'insufficient', 'YP_025292_1.insufficient_hits')
        self.assertTrue(os.path.exists(expected_output_file))

        with open(expected_output_file, 'r') as out_fh:
            seqs = list(SeqIO.parse(out_fh, 'fasta'))
            print(seqs)
            self.assertEqual(len(seqs), 5)

    def tearDown(self):
        if os.path.exists('testdir'):
            shutil.rmtree('testdir')


class TestPhylogenyPipe(BaseTestCase):
    """
    Test remaining functions used in dg core class phylogeny pipe
    """
    pass

    #def test_get_seqs_insufficient(self):
    #    """
    #    Test that dendrogenous.get_seqs() correctly renames parsed blast hit
    #    output file if there are not sufficient numbers of sequences
    #    Ensure state is correctly updated
    #    """
    #    pass
    #    #self.fail()

    #def test_align(self):
    #    """

    #    """
    #    pass
    #    #self.fail()



    #@classmethod
    #def tearDownClass(cls):
    #    """
    #    Clean up the outputs generated during testing
    #    """
    #    if os.path.exists(cls.settings['output_dir']):
    #        shutil.rmtree(cls.settings['output_dir'])

if __name__ == '__main__':
    unittest.main()
