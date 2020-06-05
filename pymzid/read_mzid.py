from time import time
from xml.etree import ElementTree
import pandas as pd
import os.path
import logging
import sys
from tqdm import tqdm

class Mzid(object):
    """
    This class holds the loaded mzIdentml object and extracts peptides and protein labels.
    Goal is to get dataframe where identified peptide sequence/z combos are linked to their m/z, protein and spectrum ID
    To do so, first pull 5 tables from the mzID file: 1.PSM, 2.Peptide, 3.ProteinGroup, 4.PeptideEvidence, 5.DBSequence

    1. PSM contains spectrum_ID, mz, z passThreshold; references to Peptide and PeptideEvidence through Pep_ref & PE_ref
    2. Peptide contains the actual amino acid sequence
    3. ProteinGroup contains passThreshold, references to DBSequence and PeptideEvidence through DBS_ref and PE_Ref
    4. PeptideEvidence contains isDecoy, references to Peptide and DBSequence through Pep_ref and DBS_ref
    5. DBSequence contains the Uniprot accession number

    From these five tables, output a peptide-centric summary, and a protein-centric summary
    Peptide-centric Summary combines PSM and Peptide, contains Sequence, mz, z, spectrumID
    Protein-centric Summary reads all 5 tables, should contain Uniprot in addition to peptide-centric summary

    """

    def __init__(self, path):
        """
        :param path: path of the mzid file to be loaded, e.g., "~/Desktop/example.mzid"


        """
        self.logger = logging.getLogger('mzid.mzid')

        self.path = os.path.join(path)
        self.root = self._parse_file()

        self.psm_df = pd.DataFrame() #self.read_psm()
        self.peptide_df = pd.DataFrame() #self.read_peptide()
        self.protein_df = pd.DataFrame() #self.read_protein()
        self.pe_df = pd.DataFrame() #self.read_pe()
        self.dbs_df = pd.DataFrame() #self.read_dbs()

        self.pep_summary_df = pd.DataFrame()



        #self.pro_summary_df = pd.DataFrame()
        #self.filtered_protein_df = pd.DataFrame()
        #self.filtered_pep_summary_df = pd.DataFrame()


    def _parse_file(self):
        """
        Get the mzid file xml root
        NB: Move time tracker to main loop when made into object oriented.

        :return: True
        """

        self.logger.info('Reading mzID file {0} as document object model...'.format(self.path))
        t1 = time()
        tree = ElementTree.parse(self.path)
        root = tree.getroot()
        # xmldoc = minidom.parse("small_test/percolator.target.mzid").childNodes[0]
        # xmldoc = minidom.parse("external_mzid_test/mzidentml-example.mzid").childNodes[0]
        # xmldoc = minidom.parse("external_mzid_test/BSA1_msgfplus_v2016_09_16.mzid").childNodes[0]

        t2 = time()
        self.logger.info('Done. Processing time: {0} seconds.'.format(round(t2 - t1, 2)))

        # Check if this is mzML 1.1 or above
        if root.attrib['version'] < '1.1.0':
            sys.exit('mzML version is not 1.1.0 or above.')

        else:
            return root


    def read_all_tables(self):
        """
        Wrapper to read all tables.
        :return:
        """
        self.read_psm()
        self.read_peptide()
        self.read_pe()
        self.read_protein()
        self.read_dbs()
        return True


    def read_psm(self):
        """
        Read in PSM table
        :return:
        """
        self.psm_df = self._read_psm()
        return None

    def read_peptide(self):
        """
        Read in Peptide table
        :return:
        """
        self.peptide_df = self._read_peptide()
        return None

    def read_protein(self):
        """
        Read in Protein table
        :return:
        """
        self.protein_df = self._read_protein()
        return None

    def read_pe(self):
        """
        Read in peptide evidence table
        """
        self.pe_df = self._read_pe()
        return None

    def read_dbs(self):
        """
        Read in dbs table
        :return:
        """
        self.dbs_df = self._read_dbs()
        return None

    def link_peptide_psm(self, take_psm_cvparams=False):
        """
        Links peptides to PSM table

        :param take_psm_cvparams: set this to true if there are additional PSM level parameters to be printed
        :return:
        """

        self.pep_summary_df = self._link_peptide_psm(self.peptide_df,
                                                     self.psm_df,
                                                     take_psm_cvparams=take_psm_cvparams)

    def _read_psm(self):
        """
        Table #1 : All PSM Identified

        Traverse to DataCollection > AnalysisData > SpectrumIdentificationList.
        The SpectrumIdentificationList (SIL) houses SpectrumIdentificationResult (SIR) elements.
        Each SIR is linked to a spectrum ID and one or more SpectrumIdentificationItem (SII) elements.
        It has the following attributes: id, spectrumID
        Each SII is a peptide-spectrum match and has the following attributes: id, rank, chargeState, peptide_ref,
        experimentalMassToCharge, calculatedMassToCharge, passThreshold

        SII has child nodes PeptideEvidenceRef and cvParam

        :return: psm_df     - A dataframe containing all peptide spectrum matches
        """

        if self.root is None:
            self.parse_file()

        list_of_all_psms = []
        root = self.root

        # Traverses down to DataCollection > AnalysisData > SpectrumIdentificationList
        data_collection = [element for element in root if element.tag.endswith('DataCollection')][0]
        analysis_data = [element for element in data_collection if element.tag.endswith('AnalysisData')][0]
        spec_id_list = [element for element in analysis_data if element.tag.endswith('SpectrumIdentificationList')][0]
        # Skip over FragmentationTable, if exists (e.g., from MSGF)
        spec_id_results = [element for element in spec_id_list if element.tag.endswith('SpectrumIdentificationResult')]

        # Traverses down to SIL > SIR > SII > PeptideEvidenceRef and cvParam
        for spec_id_result in tqdm(spec_id_results,
                                   total=len(spec_id_list),
                                   desc='Reading peptide spectrum matches'):

            spec_id_items = [element for element in spec_id_result if element.tag.endswith('SpectrumIdentificationItem')]

            # Traverse down to the SII entries inside each SIR
            for spec_id_item in spec_id_items:

                # Start an empty dictionary and add all SIR and SII level attributes
                psm_dict = {}

                psm_dict['sir_id'] = spec_id_result.attrib['id']
                psm_dict['spectrum_id'] = spec_id_result.attrib['spectrumID']

                # Get the Peptide Evidence Ref (need this to link to other tables later)
                try:
                    pep_evi_ref = [element for element in spec_id_item if element.tag.endswith('PeptideEvidenceRef')][0]
                    psm_dict['pe_id'] = pep_evi_ref.attrib['peptideEvidence_ref']
                except IndexError:
                    psm_dict['pe_id'] = None

                # Get all the SII-level attributes
                psm_dict['sii_id'] = spec_id_item.attrib['id']
                psm_dict['z'] = spec_id_item.attrib['chargeState']
                psm_dict['mz'] = spec_id_item.attrib['experimentalMassToCharge']
                psm_dict['calc_mz'] = spec_id_item.attrib['calculatedMassToCharge']
                psm_dict['pep_id'] = spec_id_item.attrib['peptide_ref']
                psm_dict['pass_threshold'] = spec_id_item.attrib['passThreshold']

                try:
                    psm_dict['rank'] = spec_id_item.attrib['rank']
                except KeyError:
                    psm_dict['rank'] = None

                # Put SIR ID, spectrum_ID, SII ID, z, m/z, theoretical m/z, Peptide ID, passThreshold, PE ID into list
                #psmList.extend([sir_id, spectrum_id, sii_id, z, mz, calc_mz, pep_id, pass_threshold, pe_id])

                # Get the list of all cvParams from each SII (e.g.,percolator scores)
                cvParams = [element for element in spec_id_item if element.tag.endswith('cvParam')]

                # Create a dictionary of name, value for each cvParam
                for cvParam in cvParams:

                    # Restrict to Top level cvParams only (no fragmentation stuff)
                    #if cvParam.parentNode != SpectrumIdentificationItem[j]:
                    #    continue

                    cvParamDict = {cvParam.attrib['name']: cvParam.attrib['value']}
                    # Add all the cvParam attributes into the psm Dictionary
                    psm_dict.update(cvParamDict)

                list_of_all_psms.append(psm_dict)

        # Convert into pandas dataframe
        psm_df = pd.DataFrame(list_of_all_psms)
        # psm_df.to_csv("all_psms.txt", sep='\t')
        # print(psm_df)
        #
        # Rename the columns in the PSM table
        #

        return psm_df

    def _read_peptide(self):
        """
        Table # 2: Peptide Sequences Identified

        ## Traverse down to Sequence Collection > Peptide > PeptideSequence and cvParams

        Get the values inside Peptide. This is needed because the SII (PSM) refers only to Peptide_ref
        Peptide_ref needs to be looked up to get the sequence.

        :return:
        """

        if self.root is None:
            self.parse_file()

        root = self.root

        list_of_all_peptides = []

        ## Traverse down to Sequence Collection > Peptide
        seq_collection = [element for element in root if element.tag.endswith('SequenceCollection')][0]
        peptides = [element for element in seq_collection if element.tag.endswith('Peptide')]

        for peptide in tqdm(peptides,
                            total=len(peptides),
                            desc='Reading peptide sequences'):

            pep_dict = {}
            pep_dict['pep_id'] = peptide.attrib['id']

            pep_seq = [element for element in peptide if element.tag.endswith('PeptideSequence')][0]


            pep_dict['seq'] = pep_seq.text # Might have to traverse down further

            # Get the cvParams from each Peptide (percolator scores)
            cvParams = [element for element in peptide if element.tag.endswith('cvParam')]

            for cvParam in cvParams:

                # Restrict to Top level cvParams only (no fragmentation stuff)
                #if cvParam.parentNode != Peptide[i]:
                #    continue

                cvParamDict = {cvParam.attrib['name']: cvParam.attrib['value']}
                # Add all the cvParam attributes into the psm Dictionary
                pep_dict.update(cvParamDict)

            list_of_all_peptides.append(pep_dict)

        peptide_df = pd.DataFrame(list_of_all_peptides)

        #
        # Rename the columns in the peptide table
        #

        return peptide_df

    def _read_protein(self):
        """
         Retrieve ProteinAmbiguityGroup for all identified proteins from the mzIdentML file.

         # Traverse to DataCollection > Analysis Data > ProteinDetectionList > Protein Ambiguity Group

        :return:
        """

        if self.root is None:
            self.parse_file()

        root = self.root

        list_of_all_proteins = []

        data_collection = [element for element in root if element.tag.endswith('DataCollection')][0]
        analysis_data = [element for element in data_collection if element.tag.endswith('AnalysisData')][0]



        # If there is no ProteinDetectionList in the mzid (apparently seen in some MSGF results)
        # Return an empty protein_df

        try:
            prot_detection_list = \
            [element for element in analysis_data if element.tag.endswith('ProteinDetectionList')][0]

        except IndexError:
            self.logger.warning('No Protein Detection List found.')
            return None


        # Traverse down to each ProteinAmbiguityGroup from each ProteinDetectionList - There are multiple PAG per PDL
        prot_ambi_groups = [element for element in prot_detection_list if element.tag.endswith('ProteinAmbiguityGroup')]

        for prot_ambi_grp in tqdm(prot_ambi_groups,
                            total=len(prot_ambi_groups),
                            desc='Reading protein ambiguity groups'):


            pag_id = prot_ambi_grp.attrib['id']

            # Traverse down to each Protein Detection Hypothesis from each Protein Ambiguity Group.
            # There may be >1 PDH per PAG in some files, but in the test files I have done.
            # (From Comet and Tide, but only single fraction) there is only 1 PDH per PAG.
            # It is possible that for multi-fraction runs there can be multiple PDH per PAG.
            prot_det_hypot = [element for element in prot_ambi_grp if element.tag.endswith('ProteinDetectionHypothesis')][0]

            # Grab the ProteinDetectionHypothesis ID, DBSequence Reference, and Threshold flag
            pdh_id = prot_det_hypot.attrib['id']
            dbs_id = prot_det_hypot.attrib['dBSequence_ref']
            pass_threshold = prot_det_hypot.attrib['passThreshold']

            # Get the cvParams from each ProteinDetectionHypothesis (e.g., percolator scores)
            cvParams = [element for element in prot_det_hypot if element.tag.endswith('cvParam')]

            cvParamsDict = {element.attrib['name']: element.attrib['value'] for element in prot_det_hypot if element.tag.endswith('cvParam')}


            # Traverse from ProteinDetectionHypothesis to each PeptideHypothesis to SpectrumIdentificationItemRef

            pep_hypots = [element for element in prot_det_hypot if element.tag.endswith('PeptideHypothesis')]

            for pep_hypot in pep_hypots:

                pe_id = pep_hypots[0].attrib['peptideEvidence_ref']

                # Traverse down to each SII Ref from each PeptideHypothesis
                spec_id_item_refs = [element for element in pep_hypot if element.tag.endswith('SpectrumIdentificationItemRef')]

                for spec_id_item_ref in spec_id_item_refs:
                    sii_id = spec_id_item_ref.attrib['spectrumIdentificationItem_ref']

                    # Put pag_id, pdh_id, dbs_ref, pass_threshold, pe_ref, and sii_ref into the dict of a particular PAG
                    protein_dict = {}

                    protein_dict['pag_id'] = pag_id
                    protein_dict['pdh_id'] = pdh_id
                    protein_dict['dbs_id'] = dbs_id
                    protein_dict['pass_threshold'] = pass_threshold
                    protein_dict['sii_id'] = sii_id

                    protein_dict.update(cvParamsDict)


                    list_of_all_proteins.append(protein_dict)

        protein_df = pd.DataFrame(list_of_all_proteins)

        return protein_df


    def _read_pe(self):
        """
        Table #4. PeptideEvidence
        This is read to see whether sequence isDecoy, and to link each PAG's PeptideHypothesis to a Peptide and a DBSequence

        Traverse down to SequenceCollection > PeptideEvidence elements
        Each PeptideEvidence refers to peptide, database sequence, and whether the peptide is decoy.

        :return:
        """

        if self.root is None:
            self.parse_file()

        root = self.root

        list_of_all_proteins = []

        seq_collection = [element for element in root if element.tag.endswith('SequenceCollection')][0]
        pep_evidences = [element for element in seq_collection if element.tag.endswith('PeptideEvidence')]

        list_of_all_pe = []

        for pep_ev in tqdm(pep_evidences,
                                   total=len(pep_evidences),
                                   desc='Reading peptide evidences'):

            pe_dict = {}

            pe_dict['pe_id'] = pep_ev.attrib['id']
            pe_dict['pep_id'] = pep_ev.attrib['peptide_ref']
            pe_dict['dbs_id'] = pep_ev.attrib['dBSequence_ref']
            pe_dict['start'] = pep_ev.attrib['start']
            pe_dict['end'] = pep_ev.attrib['end']
            pe_dict['is_decoy'] = pep_ev.attrib['isDecoy']

            list_of_all_pe.append(pe_dict)

        pe_df = pd.DataFrame(list_of_all_pe)

        return pe_df


    def _read_dbs(self):
        """
        Table #5 DBSequences
        This is read to get the Uniprot accession of each identified peptide

        :return:
        """

        if self.root is None:
            self.parse_file()

        root = self.root

        ## Traverse down to Sequence Collection > DBSequence
        seq_collection = [element for element in root if element.tag.endswith('SequenceCollection')][0]
        db_sequences = [element for element in seq_collection if element.tag.endswith('DBSequence')]

        list_of_all_dbseqs = []

        for db_seq in tqdm(db_sequences,
                                  total=len(db_sequences),
                                  desc='Reading protein database sequences'):

            db_dict = {}

            db_dict['dbs_id'] =db_seq.attrib['id']

            try:
                db_dict['length'] = db_seq.attrib['length']
            except KeyError:
                length = -1

            db_dict['accession'] =db_seq.attrib['accession']


            list_of_all_dbseqs.append(db_dict)

        dbs_df = pd.DataFrame(list_of_all_dbseqs)

        return dbs_df

    @staticmethod
    def _link_peptide_psm(peptide_df,
                          psm_df,
                          take_psm_cvparams=True):
        """
        Take the five Mzid data frames, and return a flat table

        It first combines PSM and Peptide into a peptide-centric summary, which contains:
        Peptide ID, sequence, ... (cvParams in peptide_df) ... , SIR_id, spectrum ID, SII_Id
        z, m/z, calculated m/z, pep_pass_threshold, pe_id, ... (optionally, cvParams in psm_df)

        The "take_psm_df_cvParams" flag is needed here
        because some workflows like Percolator has cvParams with identical names for different purposes
        in both the PSM section and the Peptide section.

        Percolator outputs "percolator:score",
        "percolator:Q value", and "percolator:PEP" for both Spectral Identification Item and for Peptide
        In instances where we only take the peptide level Q value but not PSM-level q value, hence if the
        take_psm_df_cvParams flag is set to False we will take only the standard columns [0:8]
        from the psm_df (without the cvParams-derived columns)
        in order to avoid conflicts in merging.

        We may have to create other flags (such as to only take certain standard
        columns from peptide_df and other dataframes) to account for other cvParams conflict.

        :param take_psm_df_cvParams: Whether to keep the cvParams columns from the PSM dataframe
        :return: True
        """

        # Here we assuming the identification scores are already sorted from top to bottom
        # #and take only one psm_id for each pep_id.
        psm_df = psm_df.groupby('pep_id', sort=False).first().reset_index().copy()

        if take_psm_cvparams is True:
            pep_summary_df = pd.merge(peptide_df, psm_df, how='left')

        elif take_psm_cvparams is False:
            pep_summary_df = pd.merge(peptide_df, psm_df[psm_df.columns[0:8]], how='left')

        return pep_summary_df


    def filter_peptide_summary(self, lysine_filter=0, protein_q=1e-2, peptide_q=1e-2, unique_only=False, require_protein_id=False):
        """
        The peptide-centric summary is then fitered by:
        - peptides that belong to any protein identified at a protein Q value
        - peptides below a certain peptide Q value
        - peptides containing certain number of lysines

        # 2017-04-06 Note the require_protein_id flag doesn't work for MSGF+ test files at the moment
        # because it seems the MSGF+ mzID files have no ProteinDetectioNList fields but instead
        # store the protein accessions inside <DBSequence>. Turn the flag off when doing MSGF+.

        :param lysine_filter: Lysine filter from command line argument
        :param protein_q: Protein-level Q value from command line argument
        :param peptide_q: Peptide-level Q value from command line argument
        :param unique_only: Only doing unique peptides
        :param require_protein_id: Require protein IDs (to filter out some mascot rows with no protein fields)
        :return: True
        """

        #
        # Filter peptides by Protein Q value.
        #
        try:
            self.filtered_protein_df = self.protein_df.loc[lambda x: x.percolator_Q_value.astype(float) < protein_q, :]
            #self.filtered_protein_df = self.filtered_protein_df.reset_index()
        except:
            print('No filtering by protein Q value done.')
            self.filtered_protein_df = self.protein_df


        #
        # Filter peptides by peptide-level Q-value filter
        #
        try:
            self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.percolator_Q_value.astype(float) < peptide_q, :]
            self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

        except:
            pass

        #
        # Filter peptides by optional lysine filter
        #

        # If the lysine filter is on and the peptide has no or more than one lysine, skup
        if lysine_filter == 1:
            try:
                self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.seq.apply(lambda y: y.count('K')) == 1, :]
                self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

            except:
                pass

        if lysine_filter == 2:
            try:
                self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.seq.apply(lambda y: y.count('K')) > 0, :]
                self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

            except:
               pass

        # 2019-05-08 Adding a new option to do only dilysine peptides in order to get RIA
        if lysine_filter == 3:
            try:
                self.pep_summary_df = self.pep_summary_df.loc[lambda x: x.seq.apply(lambda y: y.count('K')) == 2, :]
                self.pep_summary_df = self.pep_summary_df.reset_index(drop=True)

            except:
               pass

        # 2013-04-06 Again I forgot why we only chose the first five columns here. Resetting to all columns for now.

        if require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='inner') # NB: 2017-04-05 might need 'inner' here for riana to work
        elif not require_protein_id:
            self.filtered_pep_summary_df = pd.merge(self.pep_summary_df,
                                                    self.filtered_protein_df[self.filtered_protein_df.columns[[0, 5]]],
                                                    how='left')  # NB: 2017-04-05 might need 'inner' here for riana to work



        #
        # Get the protein Uniprot accession via PE and then via DBS
        #
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.pe_df, how='left')
        self.filtered_pep_summary_df = pd.merge(self.filtered_pep_summary_df, self.dbs_df, how='left')
        self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

        # Get only the peptides associated with one and only one proteins
        if unique_only:

            try:
                self.filtered_pep_summary_df = self.filtered_pep_summary_df.groupby('seq').filter(lambda x: (len(set(x['acc'])) == 1))
                self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

            except:
                pass

        self.filtered_pep_summary_df = self.filtered_pep_summary_df.reset_index(drop=True)

        return True