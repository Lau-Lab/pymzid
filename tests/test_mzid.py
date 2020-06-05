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

        global percolator_floc, msgf_floc

        percolator_floc = "tests/data/comet_percolator/percolator.target.mzid"
        msgf_floc = "tests/data/msgfplus/20180216_BSA.mzid"

        pass

    def tearDown(self):
        """

        :return:
        """
        pass


    def test_that_percolator_opens(self):
        """

        :return:
        """

        mzid = Mzid(percolator_floc)
        self.assertEqual(mzid.root.tag.endswith("MzIdentML"), True)

        return mzid

    def test_that_percolator_psm_reads(self):
        """

        :return:
        """

        mzid = self.test_that_percolator_opens()
        psm_df = mzid._read_psm()
        self.assertEqual(psm_df['pep_id'].empty, False)
        return psm_df

    def test_that_percolator_peptide_reads(self):
        """

        :return:
        """


        mzid = self.test_that_percolator_opens()
        peptide_df = mzid._read_peptide()
        self.assertEqual(peptide_df['pep_id'].empty, False)
        return peptide_df

    def test_that_percolator_peptide_psm_merges(self):
        """

        :return:
        """
        mzid = self.test_that_percolator_opens()
        psm_df = self.test_that_percolator_psm_reads()
        peptide_df = self.test_that_percolator_peptide_reads()

        merge_df = mzid._link_peptide_psm(peptide_df=peptide_df,
                                          psm_df=psm_df,
                                          take_psm_cvparams=False)

        self.assertEqual(merge_df['pe_id'].empty, False)
        # print(merge_df)

    def test_that_msgfplus_opens(self):
        """

        :return:
        """

        mzid = Mzid(msgf_floc)

        self.assertEqual(mzid.root.tag.endswith("MzIdentML"), True)
        return mzid
