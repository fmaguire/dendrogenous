#!/usr/bin/python3

import unittest
import subprocess
import dendrogenous.utils as utils
import os

from dendrogenous.test.base import BaseTestCase

class TestBinaryDependencies(BaseTestCase):
    """
    Ensure kalign, trimal and fasttree bundled dependencies run correctly
    BLASTP requires installed separately on system
    """
    def test_blastp(self):
        """
        Make sure BLASTP xml output matches what is expected when query is
        passed via stdin
        """
        test_db = os.path.join(self.test_resources,
                               "Nanoarchaeum_equitans_Kin4-M.fas")
        test_seed = (">41614809\n"
                     "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYLCNSRCTFCNIWKIY"
                     "KENPNLLKKELSKEEWIDFLSKNNKYLAWVSLTGGEPLLKDGAIDIINYLLN"
                     "NGKNVEITTNAIAYKYIIGNLEKIKNKKKLFVGVSLDGPKEIHDKLRGVPGN"
                     "FDNAIKLLEWLEQNNFNFGISVTISYQNAPYIMELFEFLEKRGWIDKVSFRV"
                     "ATISQYYRNEQKIVITKEHIQSIKEALLKIYKGYPRFRKDKFYFGNLLFLEG"
                     "KIDFDCVAGRNFLFIDPYGNVYPCLYRLDKKLGNIKENPNLLNYLGKFKRDC"
                     "KCWSECHYLPARRSSTKELLSYIWFRLFSNKKHL")

        blastp_binary_path = os.path.join(self.binary_path, 'blastp')
        blastp_cmd = ("{0} -db {1} -evalue 1e-5 -max_target_seqs 1 "
                      "-outfmt 5".format(blastp_binary_path,
                                         test_db))

        blastp_output = self.execute_cmd(blastp_cmd, input_str=test_seed)

        expected_output = self.parse_file(os.path.join(self.test_resources,
                                                       "expected_blastp_output.xml"))

        self.assertEqual(blastp_output, expected_output)


    def test_kalign(self):
        """
        Run kalign binary and ensure output alignment matches expected output
        """

        test_input = (">A\n"
                      "MMMMMSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYLCNSRCTFCNIWKIY\n"
                      ">B\n"
                      "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYYYYSRCTFCNIWKIY\n"
                      ">C\n"
                      "MSVLNKIKTVLTTPIRDIEIIEGRLKKKKYFLSLEISLTYLCNSRCTFCNIWKIY")

        expected_output = [">A",
                           "MMMMMSVLNKIKTVLTTPIRDIE---GRLKKGYYFLSLEISLTYLCNSRCTFCNIWKIY",
                           ">B",
                           "----MSVLNKIKTVLTTPIRDIE---GRLKKGYYFLSLEISLTYYYYSRCTFCNIWKIY",
                           ">C",
                           "----MSVLNKIKTVLTTPIRDIEIIEGRLKKKKYFLSLEISLTYLCNSRCTFCNIWKIY"]


        kalign_binary_path = os.path.join(self.binary_path, 'kalign')
        kalign_output = self.execute_cmd("{0}".format(kalign_binary_path),
                                         test_input)

        self.assertEqual(kalign_output, expected_output)


    def test_trimal(self):
        """
        Run trimal and ensure output mask matches expected output
        """
        trimal_binary_path = os.path.join(self.binary_path, 'trimal')
        test_input = os.path.join(self.test_resources, 'test_alignment.afa')
        trimal_cmd = "trimal -in {0} -nogaps".format(test_input)

        expected_output = [">A 52 bp",
                           "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYLCNSRCTFCNIWKIY",
                           ">B 52 bp",
                           "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYYYYSRCTFCNIWKIY",
                           ">C 52 bp",
                           "MSVLNKIKTVLTTPIRDIEGRLKKKKYFLSLEISLTYLCNSRCTFCNIWKIY"]

        trimal_output = self.execute_cmd(trimal_cmd)

        self.assertEqual(trimal_output, expected_output)


    def test_fasttree2(self):
        """
        Run fasttree2 and ensure output phylogeny matches expected output
        """
        fasttree_binary_path = os.path.join(self.binary_path, 'FastTree')
        test_mask_input = (">A 52 bp\n"
                           "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYLCNSRCTFCNIWKIY\n"
                           ">B 52 bp\n"
                           "MSVLNKIKTVLTTPIRDIEGRLKKGYYFLSLEISLTYYYYSRCTFCNIWKIY\n"
                           ">C 52 bp\n"
                           "MSVLNKIKTVLTTPIRDIEGRLKKKKYFLSLEISLTYLCNSRCTFCNIWKIY")

        fasttree_output = self.execute_cmd("{0}".format(fasttree_binary_path),
                                           input_str=test_mask_input)

        expected_output = ["(A:0.00055,B:0.06596,C:0.04382);"]

        self.assertEqual(fasttree_output, expected_output)


if __name__ == '__main__':

    unittest.main()

