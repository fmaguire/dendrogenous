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

        blastp_output = utils.execute_cmd(blastp_cmd, input_str=test_seed)

        expected_output = self.parse_file(os.path.join(self.test_resources,
                                                       "expected_blastp_output.xml"))

        print(blastp_output[-10:])
        print('\n\n')
        print(expected_output[-10:])
        self.assertEqual(blastp_output, expected_output)


    def test_kalign(self):
        """
        Run kalign binary and ensure output alignment matches expected output
        """
        self.fail()


    def test_trimal(self):
        """
        Run trimal and ensure output mask matches expected output
        """
        self.fail()


    def test_fasttree2(self):
        """
        Run fasttree2 and ensure output phylogeny matches expected output
        """
        self.fail()

if __name__ == '__main__':

    unittest.main()

