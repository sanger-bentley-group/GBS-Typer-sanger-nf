import argparse
import unittest
from unittest.mock import patch, call, ANY
from bin.mic_predictor import MicPredictor


class TestMicPredictor(unittest.TestCase):
    TEST_LANE = "26189_8#338"
    TEST_GBS_FULLGENES_RESULTS_FILE = \
        "test_data/" + TEST_LANE + "_SURFACE__fullgenes__GBS_Surface_Gene-DB_Final__results.txt"

    def setUp(self) -> None:

        self.__under_test = MicPredictor()

    def test_initialise(self) -> None:
        self.__under_test.initialise()
        self.assertEqual(self.__under_test.output, {})
        self.assertEqual(self.__under_test.pbp2x, 0.0)

    def test_get_arguments(self) -> None:
        actual = self.__under_test.get_arguments().parse_args(
            ['--input_file', 'input_file', '--pbp_file', 'pbp_file', '--output', 'outfile'])
        self.assertEqual(actual, argparse.Namespace(input_file='input_file', pbp_file='pbp_file', output='outfile'))

    def test_get_arguments_short_options(self) -> None:
        actual = self.__under_test.get_arguments().parse_args(['-i', 'input_file', '-p', 'pbp_file', '-o', 'outfile'])
        self.assertEqual(actual,
                         argparse.Namespace(
                             input_file='input_file', pbp_file='pbp_file', output='outfile'))

    def test_update_pbp_category_small_pbp2x(self) -> None:
        self.__under_test.pbp2x = 1
        self.__under_test.update_pbp_category()

        self.assertEqual(self.__under_test.output['PBP'],
                         '<=,0.5,U,<=,8.0,U,<=,0.12,S,<=,0.12,S,<=,0.12,S,NA,NA,NA,<=,0.25,S,<=,0.12,S,<=,0.12,S')

    def test_update_pbp_category_large_pbp2x(self) -> None:
        self.__under_test.pbp2x = 10
        self.__under_test.update_pbp_category()

        self.assertEqual(
            self.__under_test.output['PBP'],
            'Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,NA,NA,NA,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag'
        )

    def test_update_er_cl_category(self) -> None:
        self.__under_test.res_dict['EC'] = 'neg'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            self.__under_test.output['EC'],
            'neg,<=,0.25,S,<=,0.25,S,<=,2.0,S,<=,1.0,S,neg'
        )

    def test_update_er_cl_category(self) -> None:
        self.__under_test.res_dict['EC'] = '*R23S1*:*RPLD1*:foobar'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            self.__under_test.output['EC'],
            'neg,<=,0.25,S,<=,0.25,S,<=,2.0,S,<=,1.0,S,neg'
        )

    def test_update_gyra_parc_category(self) -> None:
        pass

    def test_update_tet_category(self) -> None:
        pass

    def test_update_other_category(self) -> None:
        pass

    @patch('bin.file_utils.FileUtils')
    @patch('bin.mic_predictor.MicPredictor.get_arguments')
    def test_run(self, mock_get_arguments, mock_file_utils) -> None:
        self.__under_test.res_dict['EC'] = '*R23S1*:*RPLD1*:foobar'
        self.__under_test.pbp2x = 2
        self.__under_test.run()


"""
    @patch('bin.process_surface_typer_results.update_protein_presence_absence')
    def test_derive_presence_absence(self, mock_update_protein_presence_absence):
        derive_presence_absence(
            self.TEST_GBS_FULLGENES_RESULTS_FILE, self.MIN_DEPTH, mock_update_protein_presence_absence)
        mock_update_protein_presence_absence.assert_has_calls([
            call('SRR1', 'SRR1-150', self.MIN_DEPTH, 138.384, featureCol, binFeatureCol, variantLookup),
            call('ALP23', 'ALP23-1', self.MIN_DEPTH, 166.103, featureCol, binFeatureCol, variantLookup),
            call('PI1', 'PI1-1', self.MIN_DEPTH, 137.787, featureCol, binFeatureCol, variantLookup),
            call('PI2A1', 'PI2A1-1', self.MIN_DEPTH, 126.262, featureCol, binFeatureCol, variantLookup),
            call('PI2A3', 'PI2A3-1', self.MIN_DEPTH, 100.967, featureCol, binFeatureCol, variantLookup),
            call('HVGA', 'HVGA1', self.MIN_DEPTH, 126.222, featureCol, binFeatureCol, variantLookup)])
"""

