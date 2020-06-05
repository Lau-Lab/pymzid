import sys
import argparse
from . import __version__

def get_args():
    """
    Get command line arguments into parser

    :return: args dict
    """

    parser = argparse.ArgumentParser(description='PyMzid {0} reads protein identification mzID files'.format(__version__))

    parser.add_argument('mzid', help='path to mzid file')
    parser.add_argument('-i', '--id', action='store_true',
                        help='only outputs rows associated with protein accession')

    parser.add_argument('-o', '--out', help='prefix of the output directory [default: mzid]',
                        default='mzid')

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()

    # Print help message if no arguments are given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    return args
