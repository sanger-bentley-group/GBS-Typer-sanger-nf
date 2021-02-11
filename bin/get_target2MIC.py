#!/usr/bin/env python3

import argparse
import sys
import re
from lib.file_utils import FileUtils


class MicPredictor:
    """ Minimum Inhibitory Concentration prediction class """

    __DEBUG = False

    OUTPUT_FIELD_SEP = ','
    RES_TARGET_DELIM = ":"
    PBP2X_DEFAULT = 'NA'
    OUTPUT_FIELDS = ['PBP', 'TET', 'EC', 'FQ', 'OTHER']

    def __init__(self):
        """ Default constructor """
        self.__res_dict = {}
        self.__output_dict = {}
        self.__pbp2X = self.PBP2X_DEFAULT
        self.actions = {
            self.update_pbp_category,
            self.update_er_cl_category,
            self.update_gyra_parc_category,
            self.update_tet_category,
            self.update_other_category
        }

    @property
    def output(self):
        return self.__output_dict

    @property
    def pbp2x(self):
        return self.__pbp2X

    @pbp2x.setter
    def pbp2x(self, value):
        self.__pbp2X = value

    @property
    def res_dict(self):
        return self.__res_dict

    @res_dict.setter
    def res_dict(self, in_dict):
        self.__res_dict = in_dict

    @staticmethod
    def search_list(arr, regex):
        """Search a text array, matching against a given regular expression"""
        for entry in arr:
            if re.search(regex, entry, re.IGNORECASE):
                return True
        return False

    def initialise(self):
        """Initialise state"""
        self.__output_dict = {}

        # Ensure ordering of res_dict as per OUTPUT_FIELDS
        for field in self.OUTPUT_FIELDS:
            self.__res_dict[field] = 'neg'

    def log(self, text):
        """Log to standard out if in debug mode"""
        if self.__DEBUG:
            print(text)

    def update_pbp_category(self):
        """Create a PBP entry in the output_dict"""
        drug = {
            'ZOX': 'Flag,Flag,Flag',
            'FOX': 'Flag,Flag,Flag',
            'TAX': 'Flag,Flag,Flag',
            'CFT': 'Flag,Flag,Flag',
            'CPT': 'Flag,Flag,Flag',
            'CZL': 'NA,NA,NA',
            'AMP': 'Flag,Flag,Flag',
            'PEN': 'Flag,Flag,Flag',
            'MER': 'Flag,Flag,Flag',
        }

        if (isinstance(self.__pbp2X, int) or isinstance(self.__pbp2X, float)) and self.__pbp2X <= 5:
            drug['ZOX'] = '<=,0.5,U'
            drug['FOX'] = '<=,8.0,U'
            drug['TAX'] = '<=,0.12,S'
            drug['CFT'] = '<=,0.12,S'
            drug['CPT'] = '<=,0.12,S'
            drug['CZL'] = 'NA,NA,NA'
            drug['AMP'] = '<=,0.25,S'
            drug['PEN'] = '<=,0.12,S'
            drug['MER'] = '<=,0.12,S'

        self.__output_dict['PBP'] = self.OUTPUT_FIELD_SEP.join(drug.values())
        self.log('PBP: ' + self.__output_dict['PBP'])

    def update_er_cl_category(self):
        """Create an EC entry in the output_dict"""
        drug = {
            'ERY': '<=,0.25,S',
            'CLI': '<=,0.25,S',
            'LZO': '<=,2.0,S',
            'SYN': '<=,1.0,S',
            'ERY_CLI': 'neg',
        }

        if self.__res_dict['EC'] != 'neg':
            res_targets = self.__res_dict['EC'].split(self.RES_TARGET_DELIM)

            if self.search_list(res_targets, r"R23S1|RPLD1|RPLV"):
                # flag everything
                drug['ERY'] = 'Flag,Flag,Flag'
                drug['CLI'] = 'Flag,Flag,Flag'
                drug['LZO'] = 'Flag,Flag,Flag'
                drug['SYN'] = 'Flag,Flag,Flag'
                drug['ERY_CLI'] = 'Flag'

            # Check for MEF
            found_mef = False
            if self.search_list(res_targets, r"MEF"):
                self.log('Found MEF')
                drug['ERY'] = '>=,1,R'
                found_mef = True

            # Check for LSA
            found_lsa = False
            if self.search_list(res_targets, r"LSA"):
                self.log('Found LSA')
                drug['CLI'] = '>=,1,R'
                found_lsa = True

            # Check for LNU
            found_lnu = False
            if self.search_list(res_targets, r"LNU"):
                self.log('Found LNU')
                drug['CLI'] = 'Flag,Flag,Flag'
                found_lnu = True

            # Check for ERM
            found_erm = False
            if self.search_list(res_targets, r"ERM"):
                self.log('Found ERM')
                drug['ERY'] = '>=,1,R'
                drug['CLI'] = '>=,1,R'
                drug['ERY_CLI'] = 'pos'
                found_erm = True

            # Check for ERY/CLI
            if found_mef and (found_lsa or found_lnu):
                self.log('Found both Mef and LSA/LNU')
                drug['ERY_CLI'] = 'pos'

            # Check for SYN
            if found_erm and found_lsa:
                self.log('Found both Erm and LSA')
                drug['SYN'] = 'Flag,Flag,Flag'
            elif found_erm and found_lnu:
                self.log('Found both Erm and LNU')
                drug['SYN'] = '<=,1.0,S'

        self.__output_dict['EC'] = self.OUTPUT_FIELD_SEP.join([self.__res_dict['EC']] + list(drug.values()))
        self.log('EC: ' + self.__output_dict['EC'])

    def update_gyra_parc_category(self):
        """Create an FQ entry in the output_dict"""
        drug = {
            'CIP': 'NA,NA,NA',
            'LFX': '<=,2,S',
        }

        if self.__res_dict['FQ'] != 'neg':
            res_targets = self.__res_dict['FQ'].split(self.RES_TARGET_DELIM)
            if self.search_list(res_targets, r"GYRA-S11L") and self.search_list(res_targets, r"PARC-S6[FY]"):
                self.log('Found GYRA-S11L:PARC-S6[FY]')
                drug['LFX'] = '>=,8,R'

            elif self.search_list(res_targets, r"PARC-D10[GY]|PARC-S6Y"):
                self.log('Found PARC-D10[GY] or PARC-S6Y')
                drug['LFX'] = '=,4,I'

            elif self.search_list(res_targets, r"PARC-D10[AN]") or self.search_list(res_targets, r"PARC-(?:D5N|S6F|S7P)"):
                self.log('Found PARC-D10[AN] or PARC-[D5N|S6F|S7P]')
                drug['LFX'] = '<=,2,S'

            else:
                drug['LFX'] = 'Flag,Flag,Flag'

        self.__output_dict['FQ'] = self.OUTPUT_FIELD_SEP.join([self.__res_dict['FQ']] + list(drug.values()))
        self.log('FQ: ' + self.__output_dict['FQ'])

    def update_tet_category(self):
        """Create a TET entry in the output_dict"""
        drug = {
            'TET': '<=,2.0,S',
        }

        if self.__res_dict['TET'] != 'neg':
            res_targets = self.__res_dict['TET'].split(self.RES_TARGET_DELIM)
            if self.search_list(res_targets, r"TET"):
                drug['TET'] = '>=,8,R'
            else:
                # New bit...
                drug['TET'] = 'neg'

        self.__output_dict['TET'] = self.OUTPUT_FIELD_SEP.join([self.__res_dict['TET']] + list(drug.values()))
        self.log('TET: ' + self.__output_dict['TET'])

    def update_other_category(self):
        """Create an OTHER entry in the output_dict"""
        drug = {
            'DAP': '<=,1,S',
            'VAN': '<=,1,S',
            'RIF': '<=,1,U',
            'CHL': '<=,4,S',
            'COT': '<=,0.5,U',
        }

        if self.__res_dict['OTHER'] != 'neg':
            res_targets = self.__res_dict['OTHER'].split(self.RES_TARGET_DELIM)
            is_new = False

            if self.__res_dict['OTHER'] != "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly" and \
                self.__res_dict['OTHER'] != "aph(3')-III:Aph3-III_AGly:Sat4A_AGly" and \
                    self.__res_dict['OTHER'] != "msr(D):MsrD_MLS":

                if not self.search_list(res_targets, r"CAT|FOLA|FOLP|RPOB|VAN"):
                    self.log('Found an ARGANNOT/RESFINDER target - flag everything')
                    # NOTE: Original perl code set SXT instead of COT, which did not get output...
                    drug['DAP'] = drug['VAN'] = drug['RIF'] = drug['CHL'] = drug['COT'] = 'Flag,Flag,Flag'
                    is_new = True

            if not is_new:
                if self.search_list(res_targets, r"CAT"):
                    drug['CHL'] = '>=,16,R'
                if self.search_list(res_targets, r"FOLA|FOLP"):
                    drug['COT'] = 'Flag,Flag,Flag'
                if self.search_list(res_targets, r"VAN"):
                    drug['VAN'] = '>=,2,U'
                if self.search_list(res_targets, r"RPOB"):
                    drug['RIF'] = 'Flag,Flag,Flag'

        self.__output_dict['OTHER'] = self.OUTPUT_FIELD_SEP.join([self.__res_dict['OTHER']] + list(drug.values()))
        self.log('OTHER: ' + self.__output_dict['OTHER'])

    def run_prediction(self, input_res_dict, pbp2X):
        """Run a MIC prediction"""
        self.initialise()

        self.__pbp2X = pbp2X

        # Copy in res dict values
        for key in input_res_dict:
            self.__res_dict[key] = input_res_dict[key]

        # Run all MIC predictor actions
        for action in self.actions:
            action()


def check_lists_equal(list_1: list, list_2: list):
    """Check if two lists contain the same items"""
    if len(list_1) != len(list_2):
        return False
    return sorted(list_1) == sorted(list_2)


def extract_pb2x_id(pbp_id_field: str):
    """Extract a pb2x_id from the given pbp field"""
    pb2x_id = MicPredictor.PBP2X_DEFAULT
    fields = pbp_id_field.split(":")
    if fields[0] != pbp_id_field:
        try:
            pb2x_id = fields[2]
        except IndexError:
            pass

    return pb2x_id


def get_arguments():
    """Parse MIC predictor command line arguments"""
    parser = argparse.ArgumentParser(description='Process MIC predictions.')
    parser.add_argument('--res_file', '-r', dest='res_file', required=True,
                        help='Input drug category allele variants file.')
    parser.add_argument('--pbp_file', '-p', dest='pbp_file', required=True,
                        help='Input pbp alleles file.')
    parser.add_argument('--output', '-o', dest='output_file', required=True,
                        help='Output filename.')
    return parser


def run():
    """Main method"""
    file_delim = '\t'
    input_res_fields = ['EC', 'FQ', 'OTHER', 'TET']
    input_pbp_fields = ['Contig', 'PBP_allele']
    mic_predictor = MicPredictor()

    # Get arguments
    args = get_arguments().parse_args()

    # Read the resistance input file
    res_headers, res_rows = FileUtils.read_delimited_id_file_with_hdrs(
        file_delim, args.res_file, len(input_res_fields)+1)
    if not check_lists_equal(res_headers, input_res_fields):
        raise RuntimeError('ERROR: File {} has an invalid format'.format(args.res_file))

    # Read the pbp input file
    pbp_headers, pbp_rows = FileUtils.read_delimited_id_file_with_hdrs(
        file_delim, args.pbp_file, len(input_pbp_fields)+1)
    if not check_lists_equal(pbp_headers, input_pbp_fields):
        raise RuntimeError('ERROR: File {} has an invalid format'.format(args.pbp_file))

    output = FileUtils.create_header_line(MicPredictor.OUTPUT_FIELDS)

    for i, sample_id in enumerate(res_rows.keys()):
        try:
            """ TODO This pb2x_id determination needs checking """
            pb2x_id = extract_pb2x_id(pbp_rows[sample_id]['PBP_allele'])
        except KeyError:
            pb2x_id = MicPredictor.PBP2X_DEFAULT
            print('WARNING: Cannot find pb2x id for sample {}'.format(sample_id))

        mic_predictor.run_prediction(res_rows[sample_id], pb2x_id)

        output += sample_id
        for col in MicPredictor.OUTPUT_FIELDS:
            output += file_delim + mic_predictor.output[col]

        if i < len(res_rows)-1:
                output += '\n'

    # Write MIC predictions
    FileUtils.write_output(output, args.output_file)


if __name__ == "__main__":
    sys.exit(run())
