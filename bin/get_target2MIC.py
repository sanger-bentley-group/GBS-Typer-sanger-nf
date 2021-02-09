#!/usr/bin/env python3

import argparse
import sys
import re
from lib.file_utils import FileUtils


class MicPredictor:
    """ Minimum Inhibitory Concentration prediction process class """

    __DEBUG = True

    INPUT_FILE_DELIM = '\t'
    OUTPUT_FIELD_SEP = ','
    RES_TARGET_DELIM = ":"
    MIC_PREDICTIONS_OUTPUT_FILE = "_mic_predictions.txt"
    NUM_RES_FILE_FIELDS = 4

    def __init__(self):
        """ Default constructor """
        self.initialise()

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
        for entry in arr:
            if re.search(regex, entry, re.IGNORECASE):
                return True
        return False

    def initialise(self):
        """Initialise variables"""
        self.__pbp2X = 0.0
        self.__output_dict = {}
        # Some new initialisation for the res dict (avoids key errors)..
        self.__res_dict = {
            'EC': 'neg',
            'FQ': 'neg',
            'OTHER': 'neg',
            'TET': 'neg',
        }

    def log(self, text):
        """Log to standard out if in debug mode"""
        if self.__DEBUG:
            print(text)

    @staticmethod
    def check_lists_equal(list_1: list, list_2: list):
        if len(list_1) != len(list_2):
            return False
        return sorted(list_1) == sorted(list_2)

    @staticmethod
    def get_arguments():
        """Parse MIC predictor command line arguments"""
        parser = argparse.ArgumentParser(description='Process MIC predictions.')
        parser.add_argument('--res_file', '-r', dest='res_file', required=True,
                            help='Input drug category allele variants file.')
        parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                            help='Output prefix for filename.')
        return parser

    def read_res_alleles_file(self, delimiter, input_filename):
        """Read the res alleles file into res_dict"""
        try:
            with open(input_filename, 'r') as fd:
                lines = fd.readlines()
                line_num = 0
                for line in lines:
                    line_num += 1

                    fields = line.strip().split(delimiter)
                    self.log("read_res_alleles_file: Read input file line: " + line)

                    if line_num == 1:
                        # read column headers
                        columns = fields

                        if not self.check_lists_equal(columns, list(self.__res_dict.keys())):
                            raise ValueError('ERROR: Input file {} has unexpected fields'.format(input_filename))
                    else:
                        field_idx = 0
                        for field in fields:
                            self.__res_dict[columns[field_idx]] = field
                            field_idx += 1

        except IOError:
            print('Cannot open filename starting "{}"'.format(input_filename))
            raise

    def update_pbp_category(self):
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

        self.__output_dict['OTHER'] = self.OUTPUT_FIELD_SEP.join([
            self.__res_dict['OTHER'],
            drug['DAP'],
            drug['VAN'],
            drug['RIF'],
            drug['CHL'],
            drug['COT']])
        self.log('OTHER: ' + self.__output_dict['OTHER'])

    def run(self):

        args = self.get_arguments().parse_args()

        self.read_res_alleles_file(self.INPUT_FILE_DELIM, args.res_file)

        # Run all MIC predictor actions
        for action in self.actions:
            action()

        output_contents = FileUtils.create_output_contents(self.output)

        # Write mic predictions
        FileUtils.write_output(output_contents, args.output + self.MIC_PREDICTIONS_OUTPUT_FILE)


if __name__ == "__main__":
    mic_predictor = MicPredictor()
    sys.exit(mic_predictor.run())
