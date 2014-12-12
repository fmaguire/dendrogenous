#!/usr/bin/env python3

import dendrogenous.core as dg
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class TestCore(unittest.TestCase):
    """
    Unit tests for the core dendrogeous class and its methods
    """
    def setUp(self):
        self.test_record = SeqRecord(\
                   Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF",
                   IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein, small")

    def test_init(self):
        """
        Ensure class init works correctly
        """
        test = dg.dendrogenous(self.test_record)

        expected_seed = ">YP_0252921_1\n
                         MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEGGEEEEVAVF"

        self.assertEqual(test.seed, expected_seed)
        self.assert(test.state, False)


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
