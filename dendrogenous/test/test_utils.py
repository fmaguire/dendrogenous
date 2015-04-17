import dendrogenous as dg
import dendrogenous.utils

from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from dendrogenous.test.base import BaseTestCase

class TestReformatAccession(BaseTestCase):

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

        truncated = dg.utils.reformat_accession(too_long)

        self.assertEqual(len(truncated.id), 20)
        self.assertEqual(truncated.id, "01234567890123456789")

    def test_reformat_accession_method_for_problematic_characters(self):
        """
        Test reformat accession works as expected
        """
        bad_char = SeqRecord(\
                   Seq("X",
                   IUPAC.protein),
                   id="|blah,|t:()",
                   name="foo",
                   description="bar, baz")

        fixed_chars = dg.utils.reformat_accession(bad_char)
        self.assertEqual(fixed_chars.id, "_blah__t___")


