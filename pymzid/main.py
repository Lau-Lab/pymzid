"""

pymzid - python mzIdentML Parser v.0.3.0.
pymzid reads in mzid files and creates flat summary tables.
Written by Edward Lau (lau1@stanford.edu) 2016-2018

Example:
    parse_mzid.py percolator.target.mzid --out=mzid.txt

"""

import os.path
import sys
import datetime
import logging

import pymzid.parse_args
from pymzid.read_mzid import Mzid
from . import __version__


def main():
    """
    Main entry point
    :return:
    """

    # Parse all the arguments
    args = pymzid.parse_args.get_args()

    # Handle command line arguments
    mzid_loc = args.mzid
    out_loc = args.out

    # Get timestamp for out files
    now = datetime.datetime.now()

    directory_to_write = os.path.join(args.out, 'mzid_' + now.strftime('%Y%m%d%H%M%S'))
    os.makedirs(directory_to_write, exist_ok=True)

    main_log = logging.getLogger('mzid')
    main_log.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(os.path.join(directory_to_write, 'mzid.log'))
    fh.setLevel(logging.DEBUG)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # add the handlers to the logger
    main_log.addHandler(fh)
    main_log.addHandler(ch)

    main_log.info(args)
    main_log.info(__version__)

    try:
        mzid = Mzid(mzid_loc)
        mzid.read_all_tables()
        mzid.link_peptide_psm()

        mzid.psm_df.to_csv(os.path.join(directory_to_write, 'mzid_psm.txt'), sep="\t")
        mzid.pe_df.to_csv(os.path.join(directory_to_write, 'mzid_pe.txt'), sep="\t")
        mzid.peptide_df.to_csv(os.path.join(directory_to_write, 'mzid_pep.txt'), sep="\t")
        if mzid.protein_df is not None:
            mzid.protein_df.to_csv(os.path.join(directory_to_write, 'mzid_prot.txt'), sep="\t")
        mzid.dbs_df.to_csv(os.path.join(directory_to_write, 'mzid_dbs.txt'), sep="\t")

        mzid.pep_summary_df.to_csv(os.path.join(directory_to_write, 'mzid_flat.txt'), sep="\t")

    except OSError as e:
        sys.exit('Failed to load mzid file. ' + str(e.errno))

    main_log.info("Success. Files saved to directory {0}".format(args.out))

    return sys.exit(os.EX_OK)

