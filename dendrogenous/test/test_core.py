#!/usr/bin/env python
from __future__ import unicode_literals
from __future__ import print_function

import dendrogenous as dg
import dendrogenous.core



import unittest
try:
    import unittest.mock as mock
except:
    import mock
import os
import shutil
import sys
import pytest
import pickle
import time

from socket import gethostname

from Bio import SeqIO
from Bio import Phylo
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
                   id="YP_025292_1", name="HokC",
                   description="toxic membrane protein, small")

        dir_paths = {"run_data": "0.run_data",
                     "blast_hits": "1.blast_hits",
                     "blast_fail": os.path.join("1.blast_hits", "insufficient_hits"),
                     "alignment": "2.alignment",
                     "mask": "3.mask",
                     "mask_fail": os.path.join("3.mask", "insufficient_sites"),
                     "tree": "4.phylogeny",
                     "named": "5.name"}

        for key in dir_paths.keys():
            dir_paths[key] = os.path.join('testdir', dir_paths[key])

        self.dir_paths = dir_paths

    def init_core():
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.logger_name = "test"
        test_class = dg.core.Dendrogenous(test_record,
                                          mock_settings)
        return test_class

    def test_init_clean(self):
        """
        Ensure class init works when provided a seqrecord and settings file
        """
        settings = mock.Mock(dg.settings.Settings)
        settings.dir_paths = self.dir_paths
        settings.logger_name = "test"

        test_class = dg.core.Dendrogenous(self.test_record,
                                          settings)

        expected_name = "YP_025292_1"
        expected_seed = (">YP_025292_1\n"
                         "MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF\n")

        self.assertEqual(test_class.seed, expected_seed)
        self.assertEqual(test_class.seq_name, expected_name)
        self.assertEqual(test_class.settings.dir_paths, self.dir_paths)

        self.assertEqual(test_class.aligned_seqs, os.path.join(self.dir_paths['alignment'],
                                                               expected_name + ".afa"))
        self.assertEqual(test_class.masked_seqs, os.path.join(self.dir_paths['mask'],
                                                              expected_name + ".mask"))
        self.assertEqual(test_class.phylogeny, os.path.join(self.dir_paths['tree'],
                                                           expected_name + ".tre"))
        self.assertEqual(test_class.named_phylogeny, os.path.join(self.dir_paths['named'],
                                                                  expected_name + ".named_tre"))

    def test_init_bad_seqrecord(self):
        """
        Ensure ValueError is raised if class is instantiated without valid seqrecord
        """
        settings = mock.Mock(dg.settings.Settings)

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
                   id="YP_025292_1", name="HokC",
                   description="toxic membrane protein, small")

        dir_paths = {"run_data": "0.run_data",
                     "blast_hits": "1.blast_hits",
                     "blast_fail": os.path.join("1.blast_hits", "insufficient_hits"),
                     "alignment": "2.alignment",
                     "mask": "3.mask",
                     "mask_fail": os.path.join("3.mask", "insufficient_sites"),
                     "tree": "4.phylogeny",
                     "named": "5.name"}

        for key in dir_paths.keys():
            dir_paths[key] = os.path.join('testdir', dir_paths[key])

        self.dir_paths = dir_paths


    def test__blast_runs(self):
        """
        Make sure the __blast method correctly runs and returns decoded xml output
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.blast_settings = {'num_seqs': 1,
                                        'evalue': 1e-5}
        mock_settings.logger_name = "test"
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
        mock_settings.dir_paths = self.dir_paths
        mock_settings.logger_name = "test"

        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.dbconfig = pickle.load(secret)

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
        mock_settings.dir_paths = self.dir_paths
        mock_settings.logger_name = "test"

        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.dbconfig = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        blastp_xml = self.parse_file(os.path.join(self.test_resources,
                                                  "broken_blastp_output.xml"))
        blastp_output = "\n".join(blastp_xml)

        parsed_blast = test_class._parse_blast(blastp_output)
        expected_output = []

        self.assertEqual(parsed_blast, expected_output)


    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_get_seqs(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory - integration test with _blast and _parse
        """

        #configure all the dependencies in a mock settings object
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.minimums = {'min_seqs': 3}
        mock_settings.blast_settings = {'num_seqs': 2,
                                        'evalue': 5}
        mock_settings.genomes= [os.path.join(self.test_resources, "Escherichia_coli_O157_H7_str._Sakai.fas"),
                                os.path.join(self.test_resources, "Escherichia_coli_IAI39.fas"),
                                os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")]
        mock_settings.logger_name = "test"

        # make output dir that is normally done by runner
        os.mkdir('testdir')
        os.mkdir(mock_settings.dir_paths['blast_hits'])

        # load db settings from secret pickle file
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.dbconfig = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        test_class.get_seqs()

        expected_output_file = os.path.join(self.dir_paths['blast_hits'], 'YP_025292_1.fas')
        self.assertTrue(os.path.exists(expected_output_file))

        with open(expected_output_file, 'r') as out_fh:
            seqs = list(SeqIO.parse(out_fh, 'fasta'))
            print(seqs)
            self.assertEqual(len(seqs), 6)


    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_get_seqs_fails_correctly(self):
        """
        Test dendrogenous.get_seqs() works correctly and writes a file of seqs
        to the appropriate directory - integration test with _blast and _parse
        Test it correctly identifies too few hits and moves file to insufficient hits dir
        """

        #configure all the dependencies in a mock settings object
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'blastp': os.path.join(self.binary_path, "blastp")}
        mock_settings.minimums = {'min_seqs': 10}
        mock_settings.blast_settings = {'num_seqs': 2,
                                        'evalue': 5}
        mock_settings.genomes= [os.path.join(self.test_resources, "Escherichia_coli_O157_H7_str._Sakai.fas"),
                                os.path.join(self.test_resources, "Escherichia_coli_IAI39.fas"),
                                os.path.join(self.test_resources, "Nanoarchaeum_equitans_Kin4-M.fas")]
        mock_settings.logger_name = "test"

        # make output dir that is normally done by runner
        os.mkdir('testdir')
        os.mkdir(mock_settings.dir_paths['blast_hits'])
        os.mkdir(mock_settings.dir_paths['blast_fail'])

        # load db settings from secret pickle file
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.dbconfig = pickle.load(secret)

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with self.assertRaises(dg.utils.GetSeqFail):
            test_class.get_seqs()

        expected_output_file = os.path.join(self.dir_paths['blast_fail'], 'YP_025292_1.insufficient_hits')
        self.assertTrue(os.path.exists(expected_output_file))

        with open(expected_output_file, 'r') as out_fh:
            seqs = list(SeqIO.parse(out_fh, 'fasta'))
            print(seqs)
            self.assertEqual(len(seqs), 6)

    def tearDown(self):
        if os.path.exists('testdir'):
            shutil.rmtree('testdir')


class TestPhylogenyPipe(BaseTestCase):
    """
    Test remaining functions used in dg core class phylogeny pipe
    """
    def setUp(self):
        self.test_record = SeqRecord(\
                   Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF",
                   IUPAC.protein),
                   id="test")
        self.test_dir = 'testdir2'

        dir_paths = {"run_data": "0.run_data",
                     "blast_hits": "1.blast_hits",
                     "blast_fail": os.path.join("1.blast_hits", "insufficient_hits"),
                     "alignment": "2.alignment",
                     "mask": "3.mask",
                     "mask_fail": os.path.join("3.mask", "insufficient_sites"),
                     "tree": "4.phylogeny",
                     "named": "5.name"}

        for key in dir_paths.keys():
            dir_paths[key] = os.path.join(self.test_dir, dir_paths[key])

        self.dir_paths = dir_paths
        os.mkdir(self.test_dir)

    def test_align(self):
        """
        Check align runs correctly
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'kalign': os.path.join(self.binary_path, "kalign")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['blast_hits'])
        os.mkdir(self.dir_paths['alignment'])

        shutil.copy(os.path.join(self.test_resources, 'test.fas'), self.dir_paths['blast_hits'])

        expected_file = os.path.join(self.dir_paths['alignment'], 'test.afa')
        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with mock.patch.object(test_class, 'get_seqs', return_value=None) as mock_method:
            test_class.align()

        self.assertEqual(self.file_hash(expected_file),
                         self.file_hash(os.path.join(self.test_resources, 'test_alignment.afa')))
        self.assertFalse(mock_method.called)


    def test_align_calls_seqs_if_seqs_missing(self):
        """
        Check align runs called get_seqs if seqs file is missing
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'kalign': os.path.join(self.binary_path, "kalign")}
        mock_settings.logger_name = "test"


        os.mkdir(self.dir_paths['blast_hits'])
        os.mkdir(self.dir_paths['alignment'])

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        # patching out get_seqs so this doesn't run and the output check function as this will
        # also fail due to get_seqs not actually running and thus outputting anything
        with mock.patch.object(test_class, 'get_seqs', return_value=None) as mock_method:
            with mock.patch.object(test_class, '_check_output', return_value=None) as mock_check:
                test_class.align()

        self.assertTrue(mock_method.called)
        self.assertTrue(mock_check.called)

    def test_mask_normal(self):
        """
        Check mask runs correctly when alignment file exists and mask has enough sites
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.minimums = {'min_sites': 29}
        mock_settings.binary_paths = {'trimal': os.path.join(self.binary_path, "trimal")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['alignment'])
        os.mkdir(self.dir_paths['mask'])

        shutil.copy(os.path.join(self.test_resources, 'test_alignment.afa'),
                    os.path.join(self.dir_paths['alignment'], 'test.afa'))

        expected_file = os.path.join(self.dir_paths['mask'], 'test.mask')
        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with mock.patch.object(test_class, 'align', return_value=None) as mock_method:
            test_class.mask()

        self.assertEqual(self.file_hash(expected_file),
                         self.file_hash(os.path.join(self.test_resources, 'test_mask.mask')))
        self.assertFalse(mock_method.called)


    def test_mask_needs_automated_mask(self):
        """
        Check mask correctly reruns trimal with automated if nogaps produces too small a mask
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.minimums = {'min_sites': 40}
        mock_settings.binary_paths = {'trimal': os.path.join(self.binary_path, "trimal")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['alignment'])
        os.mkdir(self.dir_paths['mask'])

        shutil.copy(os.path.join(self.test_resources, 'test_alignment_auto_mask.afa'),
                    os.path.join(self.dir_paths['alignment'], 'test.afa'))

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        expected_file = os.path.join(self.dir_paths['mask'], 'test.mask')


        with mock.patch.object(test_class, 'align', return_value=None) as mock_method:
            test_class.mask()

        self.assertTrue(os.path.exists(expected_file))
        self.assertEqual(self.file_hash(expected_file),
                         self.file_hash(os.path.join(self.test_resources, 'test_mask_automated.mask')))
        self.assertFalse(mock_method.called)


    def test_mask_fails_correctly(self):
        """
        Ensure mask fails correctly if automated and nogaps masking
        still results in too short a file
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.logger_name = "test"
        mock_settings.minimums = {'min_sites': 100}
        mock_settings.binary_paths = {'trimal': os.path.join(self.binary_path, "trimal")}

        os.mkdir(self.dir_paths['alignment'])
        os.mkdir(self.dir_paths['mask'])
        os.mkdir(self.dir_paths['mask_fail'])

        shutil.copy(os.path.join(self.test_resources, 'test_alignment.afa'),
                    os.path.join(self.dir_paths['alignment'], 'test.afa'))

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with self.assertRaises(dg.utils.MaskFail):
            test_class.mask()

        not_expected_file = os.path.join(self.dir_paths['mask'], 'test.mask')
        expected_file = os.path.join(self.dir_paths['mask_fail'], 'test.mask_too_short')
        self.assertFalse(os.path.exists(not_expected_file))
        self.assertTrue(os.path.exists(expected_file))

    def test_mask_calls_align_if_alignment_missing(self):
        """
        Ensure dg.mask() calls dg.align() if it can't find the alignment file
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.minimums = {'min_sites': 29}
        mock_settings.binary_paths = {'trimal': os.path.join(self.binary_path, "trimal")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['alignment'])
        os.mkdir(self.dir_paths['mask'])

        expected_file = os.path.join(self.dir_paths['mask'], 'test.mask')
        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        # mock align but mask will fail and align is not actually aligning thus
        # the try/except hack
        with mock.patch.object(test_class, 'align', return_value=None) as mock_method:
            try:
                test_class.mask()
            except:
                pass

        self.assertTrue(mock_method.called)

    # not sure how to mock masked_seqs to return a given value
    #def test__mask_check(self):
    #    """
    #    ensure mask check returns correct sequence length for first seq in masked
    #    fasta
    #    """

    #    mock_settings = mock.Mock(dg.settings.Settings)
    #    mock_settings.dir_paths = self.dir_paths
    #    test_class = dg.core.Dendrogenous(self.test_record,
    #                                      mock_settings)

    #    test_file = os.path.join(self.test_resources, 'test_mask.mask')

    #    with mock.patch.object(test_class, 'masked_seqs', new_callable=mock.PropertyMock) as mock_attr:
    #        mock_attr.return_value = test_file
    #        mask_size = test_class._mask_check()

    #    self.assertEqual(mask_size, 52)

    def test_phylogeny_normal(self):
        """
        Ensure phylogeny code correctly runs and generates of phylogeny in the right place
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'FastTree': os.path.join(self.binary_path, "FastTree")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['mask'])
        os.mkdir(self.dir_paths['tree'])

        shutil.copy(os.path.join(self.test_resources, 'test_mask.mask'),
                    os.path.join(self.dir_paths['mask'], 'test.mask'))

        expected_file = os.path.join(self.dir_paths['tree'], 'test.tre')
        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with mock.patch.object(test_class, 'mask', return_value=None) as mock_method:
            test_class.estimate_phylogeny()

        self.assertEqual(self.file_hash(expected_file),
                         self.file_hash(os.path.join(self.test_resources, 'test_tree.tre')))
        self.assertFalse(mock_method.called)


    def test_phylogeny_runs_mask_if_mask_is_missing(self):
        """
        Ensure phylogeny code correctly calls self.mask() if phylogeny is missing
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.binary_paths = {'FastTree': os.path.join(self.binary_path, "FastTree")}
        mock_settings.logger_name = "test"

        os.mkdir(self.dir_paths['mask'])
        os.mkdir(self.dir_paths['tree'])

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        with mock.patch.object(test_class, 'mask', return_value=None) as mock_method:
            with mock.patch.object(test_class, '_check_output', return_value=None) as mock_check:
                test_class.estimate_phylogeny()

        self.assertTrue(mock_method.called)
        self.assertTrue(mock_check.called)


    @pytest.mark.skipif("gethostname() != 'zorya'")
    def test_phylogeny_rename(self):
        """
        Ensure phylogeny rename works as expected
        """
        mock_settings = mock.Mock(dg.settings.Settings)
        mock_settings.dir_paths = self.dir_paths
        mock_settings.logger_name = "test"
        with open(os.path.join(self.test_resources, 'secret.pickle'), 'rb') as secret:
            mock_settings.dbconfig = pickle.load(secret)

        os.mkdir(self.dir_paths['tree'])
        os.mkdir(self.dir_paths['named'])

        shutil.copy(os.path.join(self.test_resources, 'name_test.tre'),
                    os.path.join(self.dir_paths['tree'], 'test.tre'))

        test_class = dg.core.Dendrogenous(self.test_record,
                                          mock_settings)

        expected_file = os.path.join(self.dir_paths['named'], 'test.named_tre')
        with mock.patch.object(test_class, 'estimate_phylogeny', return_value=None) as mock_method:
            test_class.name_phylogeny()

        self.assertFalse(mock_method.called)
        self.assertTrue(os.path.exists(expected_file))

        output_tree = Phylo.read(expected_file, 'newick')
        # this parser sucks it only returns the last whitespace separated label
        terminal_labels = [leaf.name for leaf in output_tree.get_terminals()]

        self.assertEqual(len(terminal_labels), 10)
        self.assertIn("TAXA", terminal_labels)
        self.assertIn("SEQUENCE", terminal_labels)


    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

class TestBuildTraining(BaseTestCase):

    def test_function(self):

        training_dir = os.path.join(self.test_resources, "training")
        settings = {"class_defs": {"host": ["Alveolata",
                            "Paramecium",
                            "Tetrahymena",
                            "Oxytricha"],
                   "endosymbiont": ["Arabidopsis",
                            "Chlamydomonas",
                            "Ostreococcus",
                            "Micromonas",
                            "Chlorella",
                            "Physcomitrella"],
                   "food": ["Bacteria",
                            "Bacillus",
                            "Escherichia",
                            "Salmonella",
                            "Chlamydophila",
                            "Chlorobium",
                            "Deinococcus",
                            "Caulobacter"],
                   "unknown":["Saccharomyces",
                            "Neurospora",
                            "Homo",
                            "Mus",
                            "Dictyostelium",
                            "Toxoplasma",
                            "Guillardia",
                            "Bigelowiella",
                            "Emiliania",
                            "Aureococcus",
                            "Ectocarpus",
                            "Schizosaccharomyces",
                            "Amycolatopsis",
                            "Aquifex",
                            "Sulfolobus",
                            "Nanoarchaeum",
                            "Haloferax",
                            "Methanococcus",
                            "Cenarchaeum"]},
                "class_locs": {"host": os.path.join(training_dir, "host"),
                   "endosymbiont": os.path.join(training_dir, "endosymbiont"),
                   "food": os.path.join(training_dir, "food"),
                   "unknown": os.path.join(training_dir, "unknown")}}

        label_def = settings['class_defs']
        label_loc = settings['class_locs']

        a = dg.core.BuildTraining(label_def, label_loc)

        X, y, encoded_labels = a.build_training()


        self.assertEqual(X.shape, (46,4))

        self.assertEqual(y.shape, (46,1))

        self.assertEqual(type(encoded_labels), dict)


if __name__ == '__main__':
    unittest.main()
