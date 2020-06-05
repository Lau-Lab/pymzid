"""

pymzid - python mzIdentML Parser v.0.3.0.
pymzid reads in mzid files and creates flat summary tables.
Written by Edward Lau (lau1@stanford.edu) 2016-2018

Example:
    parse_mzid.py percolator.target.mzid --out=mzid.txt

"""

import os.path
import sys
import pymzid.parse_args
from pymzid.read_mzid import Mzid


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

    try:
        mzid = Mzid(mzid_loc)
        mzid.merge_tables(take_psm_df_cvParams=True)

        #mzid.filter_peptide_summary(lysine_filter=0, protein_q=1, peptide_q=1, unique_only=False, require_protein_id=False)
        #mzid.filtered_pep_summary_df.to_csv(out_loc, sep='\t')

        mzid.psm_df.to_csv(out_loc + 'psm.txt', sep='\t')
        mzid.pe_df.to_csv(out_loc + 'pe.txt', sep='\t')
        mzid.peptide_df.to_csv(out_loc + 'peptide.txt', sep='\t')
        if mzid.protein_df is not None:
            mzid.protein_df.to_csv(out_loc + 'proteins.txt', sep='\t')
        mzid.dbs_df.to_csv(out_loc + 'dbs.txt', sep='\t')

        mzid.pep_summary_df.to_csv(out_loc + '.txt', sep='\t')

    except OSError as e:
        sys.exit('Failed to load mzid file. ' + str(e.errno))

    return sys.exit(os.EX_OK)

