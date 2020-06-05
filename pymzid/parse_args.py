import sys
import argparse

def get_args():
    """
    Get command line arguments into parser

    :return: args dict
    """

    parser = argparse.ArgumentParser(description='PyMzid v.0.3.0 reads protein identification mzID files')

    parser.add_argument('mzid', help='path to mzid file')
    parser.add_argument('-i', '--id', action='store_true',
                        help='only outputs rows associated with protein accession')

    # TODO: this needs to be a directory and adds in time like in riana
    parser.add_argument('-o', '--out', help='name of the output files [default: mzid.txt]',
                        default='mzid')

    args = parser.parse_args()

    # Print help message if no arguments are given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    return args
