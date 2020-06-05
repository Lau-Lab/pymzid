import unittest
from pymzid.read_mzid import Mzid

class MzidTest(unittest.TestCase):
    """
    Test cases involving reading mzid files
    """

    def setUp(self):
        """

        :return:
        """

        global percolator_floc

        percolator_floc = "tests/data/comet_percolator/percolator.target.mzid"

        pass

    def tearDown(self):
        """

        :return:
        """
        pass


    def test_that_percolator_opens(self):

        mzid = Mzid(percolator_floc)
        mzid.parse_file()

        return mzid

    def test_that_percolator_psm_reads(self):

        mzid = self.test_that_percolator_opens()

        psm_df = mzid._read_psm()

        print(psm_df)

        self.assertEqual(psm_df['pep_id'] is not None)

    def test_thatpercolator_peptide_reads(self):

        mzid = self.test_that_percolator_opens()
        peptide_df = mzid._read_peptide()

        print(peptide_df)

        self.assertEqual(peptide_df['pep_id'] is not None)

