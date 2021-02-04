import argparse
import sys
import re
from bin.file_utils import FileUtils


class MicPredictor:
    """ Minimum Inhibitory Concentration prediction process class """

    __DEBUG = True

    INPUT_FILE_DELIM = '\t'
    OUTPUT_FIELD_SEP = ','
    RES_TARGET_DELIM = ":"

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
    def search_array(arr, regex):
        for entry in arr:
            if re.search(regex, entry, re.IGNORECASE):
                return True
        return False

    def initialise(self):
        """Initialise variables"""
        self.__pbp2X = 0.0
        self.__drug = {}
        self.__res_dict = {}
        self.__output_dict = {}

    def log(self, text):
        """Log to standard out if in debug mode"""
        if self.__DEBUG:
            print(text)

    @staticmethod
    def get_arguments():
        """Parse MIC predictor command line arguments"""
        parser = argparse.ArgumentParser(description='Process MIC predictions.')
        parser.add_argument('--input_file', '-i', dest='input_file', required=True,
                            help='Input res file.')
        parser.add_argument('--pbp_file', '-p', dest='pbp_file', required=True,
                            help='Input pbp file.')
        parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                            help='Output prefix for filename.')
        return parser

    def update_pbp_category(self):
        self.__drug['ZOX'] = 'Flag,Flag,Flag'
        self.__drug['FOX'] = 'Flag,Flag,Flag'
        self.__drug['TAX'] = 'Flag,Flag,Flag'
        self.__drug['CFT'] = 'Flag,Flag,Flag'
        self.__drug['CPT'] = 'Flag,Flag,Flag'
        self.__drug['CZL'] = 'NA,NA,NA'
        self.__drug['AMP'] = 'Flag,Flag,Flag'
        self.__drug['PEN'] = 'Flag,Flag,Flag'
        self.__drug['MER'] = 'Flag,Flag,Flag'

        if (isinstance(self.__pbp2X, int) or isinstance(self.__pbp2X, float)) and self.__pbp2X <= 5:
            self.__drug['ZOX'] = '<=,0.5,U'
            self.__drug['FOX'] = '<=,8.0,U'
            self.__drug['TAX'] = '<=,0.12,S'
            self.__drug['CFT'] = '<=,0.12,S'
            self.__drug['CPT'] = '<=,0.12,S'
            self.__drug['CZL'] = 'NA,NA,NA'
            self.__drug['AMP'] = '<=,0.25,S'
            self.__drug['PEN'] = '<=,0.12,S'
            self.__drug['MER'] = '<=,0.12,S'

        self.__output_dict['PBP'] = self.OUTPUT_FIELD_SEP.join(self.__drug.values())
        self.log('PBP: ' + self.__output_dict['PBP'])

    def update_er_cl_category(self):
        self.__drug['ERY'] = '<=,0.25,S'
        self.__drug['CLI'] = '<=,0.25,S'
        self.__drug['SYN'] = '<=,1.0,S'
        self.__drug['LZO'] = '<=,2.0,S'
        self.__drug['ERY_CLI'] ='neg'

        if self.__res_dict['EC'] == 'neg':
            self.__output_dict['EC'] = self.OUTPUT_FIELD_SEP.join([
                self.__res_dict['EC'],
                self.__drug['ERY'],
                self.__drug['CLI'],
                self.__drug['LZO'],
                self.__drug['SYN'],
                self.__drug['ERY_CLI']])
        else:
            res_targets = self.__res_dict['EC'].split(self.RES_TARGET_DELIM)

            if self.search_array(res_targets, r"R23S1|RPLD1|RPLV"):
                # flag everything
                self.__drug['ERY'] = 'Flag,Flag,Flag'
                self.__drug['CLI'] = 'Flag,Flag,Flag'
                self.__drug['SYN'] = 'Flag,Flag,Flag'
                self.__drug['LZO'] = 'Flag,Flag,Flag'
                self.__drug['ERY_CLI'] = 'Flag'

            # Check for MEF
            found_mef = False
            if self.search_array(res_targets, r"MEF"):
                self.log('Found MEF')
                self.__drug['ERY'] = '>=,1,R'
                found_mef = True

            # Check for LSA
            found_lsa = False
            if self.search_array(res_targets, r"LSA"):
                self.log('Found LSA')
                self.__drug['CLI'] = '>=,1,R'
                found_lsa = True

            # Check for LNU
            found_lnu = False
            if self.search_array(res_targets, r"LNU"):
                self.log('Found LNU')
                self.__drug['CLI'] = 'Flag,Flag,Flag'
                found_lnu = True

            # Check for ERM
            found_erm = False
            if self.search_array(res_targets, r"ERM"):
                self.log('Found ERM')
                self.__drug['ERY'] = '>=,1,R'
                self.__drug['CLI'] = '>=,1,R'
                self.__drug['ERY_CLI'] = 'pos'
                found_erm = True

            # Check for ERY/CLI
            if found_mef and (found_lsa or found_lnu):
                self.log('Found both Mef and LSA/LNU')
                self.__drug['ERY_CLI'] = 'pos'

            # Check for SYN
            if found_erm and found_lsa:
                self.log('Found both Erm and LSA')
                self.__drug['SYN'] = 'Flag,Flag,Flag'
            elif found_erm and found_lnu:
                self.log('Found both Erm and LNU')
                self.__drug['SYN'] = '<=,1.0,S'

            self.__output_dict['EC'] = self.OUTPUT_FIELD_SEP.join(
                [self.__res_dict['EC'], self.__drug['ERY'], self.__drug['CLI'], self.__drug['LZO'], self.__drug['SYN'], self.__drug['ERY_CLI']])

        self.log('EC: ' + self.__output_dict['EC'])

    def update_gyra_parc_category(self):
        self.__drug['LFX'] = '<=,2,S'
        self.__drug['CIP'] = 'NA,NA,NA'

        if self.__res_dict['FQ'] != 'neg':
            res_targets = self.__res_dict['FQ'].split(self.RES_TARGET_DELIM)
            if self.search_array(res_targets, r"GYRA-S11L") and self.search_array(res_targets, r"PARC-S6[FY]"):
                self.log('Found GYRA-S11L:PARC-S6[F|Y]')
                self.__drug['LFX'] = '>=,8,R'
            elif self.search_array(res_targets, r"PARC-D10[GY]|PARC-S6Y"):
                self.log('Found PARC-D10[G|Y] or PARC-S6Y')
                self.__drug['LFX'] = '=,4,I'
            elif self.search_array(res_targets, r"PARC-D10[AN]") or self.search_array(res_targets, r"PARC-[D5N|S6F|S7P]"):
                self.log('Found PARC-D10[A|N] or PARC-[D5N|S6F|S7P]')
                self.__drug['LFX'] = '<=,2,S'
            else:
                self.__drug['LFX'] = 'Flag,Flag,Flag'

        self.__output_dict['FQ'] = \
                self.OUTPUT_FIELD_SEP.join([self.__res_dict['FQ'], self.__drug['CIP'], self.__drug['LFX']])
        self.log('FQ: ' + self.__output_dict['FQ'])

    def update_tet_category(self):
        self.__drug['TET'] = '<=,2.0,S'

        if self.__res_dict['TET'] != 'neg':
            res_targets = self.__res_dict['TET'].split(":")
            if self.search_array(res_targets, r"TET"):
                self.__drug['TET'] = '>=,8,R'

        self.__output_dict['TET'] = self.OUTPUT_FIELD_SEP.join([self.__res_dict['TET'], self.__drug['TET']])
        self.log('TET: ' + self.__output_dict['TET'])

    def update_other_category(self):
        self.__drug['DAP'] = '<=,1,S'
        self.__drug['VAN'] = '<=,1,S'
        self.__drug['RIF'] = '<=,1,U'
        self.__drug['CHL'] = '<=,4,S'
        self.__drug['COT'] = '<=,0.5,U'

        if self.__res_dict['OTHER'] != 'neg':
            res_targets = self.__res_dict['OTHER'].split(":")
            is_new = False

            if self.__res_dict['OTHER'] != "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly" and \
                self.__res_dict['OTHER'] != "aph(3')-III:Aph3-III_AGly:Sat4A_AGly" and \
                self.__res_dict['OTHER'] != "msr(D):MsrD_MLS":

                if not self.search_array(res_targets, r"CAT|FOLA|FOLP|RPOB|VAN"):
                    self.log('Found an ARGANNOT/RESFINDER target - flag everything')
                    self.__drug['DAP'] = self.__drug['VAN'] = self.__drug['RIF'] = self.__drug['CHL'] = \
                    self.__drug['SXT'] = 'Flag,Flag,Flag'
                    is_new = True

            if not is_new:
                if self.search_array(res_targets, r"CAT"):
                    self.__drug['CHL'] = '>=,16,R'
                if self.search_array(res_targets, r"FOLA|FOLP"):
                    self.__drug['COT'] = 'Flag,Flag,Flag'
                if self.search_array(res_targets, r"VAN"):
                    self.__drug['VAN'] = '>=,2,U'
                if self.search_array(res_targets, r"RPOB"):
                    self.__drug['RIF'] = 'Flag,Flag,Flag'

        self.__output_dict['OTHER'] = self.OUTPUT_FIELD_SEP.join([
            self.__res_dict['OTHER'],
            self.__drug['DAP'],
            self.__drug['VAN'],
            self.__drug['RIF'],
            self.__drug['CHL'],
            self.__drug['COT']])
        self.log('OTHER: ' + self.__output_dict['OTHER'])

    def run(self):

        args = self.get_arguments().parse_args()

        self.__res_dict = FileUtils.read_delimited_keyvalue_file(
            self.INPUT_FILE_DELIM, args.input_file)

        # TODO read input pbp file.

        # Run all MIC predictor actions
        for action in self.actions:
            action()

        output_contents = FileUtils.create_output_contents(self.output)

        # Write mic predictions
        FileUtils.write_output(output_contents, args.output + "_mic_predictions.txt")


def main():
    mic_predictor = MicPredictor
    mic_predictor.run()


if __name__ == "__main__":
    sys.exit(main())
