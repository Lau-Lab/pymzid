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

        pass