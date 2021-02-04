import unittest
import argparse
from unittest.mock import patch, call, ANY

from bin.combine_results import get_content_with_id, get_sero_res_contents, get_arguments, main


class TestCombineResults(unittest.TestCase):

    TEST_LANE = '25292_2#85'
    TEST_DATA_SEROTYPE = 'test_data/input/' + TEST_LANE + '_SeroType_Results.txt'
    TEST_DATA_RES_INCIDENCE = 'test_data/input/' + TEST_LANE + '_res_incidence.txt'
    TEST_DATA_RES_ALLELES = 'test_data/input/' + TEST_LANE + '_res_alleles.txt'
    TEST_DATA_RES_VARIANTS = 'test_data/input/' + TEST_LANE + '_res_gbs_variants.txt'
    TEST_DATA_PBP_ALLELE = 'test_data/input/' + TEST_LANE + '_GBS1A-1_PBP_existing_allele.txt'
    TEST_OUTPUT = "test_data/output/" + TEST_LANE + "_output.txt"

    def test_should_get_content_with_id_alleles(self):
        actual = get_content_with_id(self.TEST_LANE, self.TEST_DATA_RES_ALLELES)
        self.assertEqual(actual, "ID\tEC\tFQ\tOTHER\tTET\n25292_2#85\tERMB(ERMB-1):MEF(MEF-1):23S1:23S3\tPARC:GYRA-V1A,M2Q,G3W,K4W\tMsrD_MLS(MsrD_296):Ant6-Ia_AGly(Ant6-Ia_1633):Sat4A_AGly(Sat4A_586):Aph3-III_AGly(Aph3-III_268):msr(D)(msr(D)_2):ant(6)-Ia(ant(6)-Ia):aph(3')-III(aph(3')-III)\tTETO(TETO-1)\n")

    def test_should_get_content_with_id_variants(self):
        actual = get_content_with_id(self.TEST_LANE, self.TEST_DATA_RES_VARIANTS)
        self.assertEqual(actual, "ID\t23S1\t23S3\tGYRA\tPARC\tRPOBGBS-1\tRPOBGBS-2\tRPOBGBS-3\tRPOBGBS-4\n25292_2#85\t23S1\t23S3\tGYRA-V1A,M2Q,G3W,K4W\tPARC\t\t\t\t\n")

    def test_should_get_sero_res_contents(self):
        actual = get_sero_res_contents(self.TEST_LANE, self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE)
        self.assertEqual(actual, 'ID\tSerotype\t23S1\t23S3\tCAT\tERM\tGYRA\tLNU\tLSA\tMEF\tPARC\tRPOBGBS-1\tRPOBGBS-2\tRPOBGBS-3\tRPOBGBS-4\tTET\n25292_2#85\tIII\t+\t+\t-\t+\t+\t-\t-\t+\t+\t-\t-\t-\t-\t+\n')

    def test_should_get_pbp_contents(self):
        actual = get_content_with_id(self.TEST_LANE, self.TEST_DATA_PBP_ALLELE)
        self.assertEqual(actual, 'ID\tContig\tPBP_allele\n25292_2#85\t.25292_2_85.9:39495-40455(+)\t2||GBS_1A\n')

    def test_arguments_sero_res(self):
        actual = get_arguments().parse_args(['sero_res', '--id', 'id', '--serotyper_results', 'sero_file', '--res_incidence_results', 'res_file',
                                            '--res_alleles_results', 'allele_file', '--res_variants_results', 'variants_file', '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='sero_res', id='id', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))

    def test_arguments_surface_typing(self):
        actual = get_arguments().parse_args(['surface_typer', '--id', 'id', '--surface_incidence_results', 'file1', '--surface_variants_results', 'file2',
                                             '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='surface_typer', id='id', surface_inc='file1', surface_variants='file2', output='output_prefix'))

    def test_arguments_pbp_typing(self):
        actual = get_arguments().parse_args(['pbp_typer', '--id', 'id', '--pbp_existing_allele_results', 'pbp_file', '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='pbp_typer', id='id', pbp_allele='pbp_file', output='output_prefix'))

    def test_arguments_short_options_sero_res(self):
        actual = get_arguments().parse_args(['sero_res', '-i', 'id', '-s', 'sero_file', '-r', 'res_file', '-a', 'allele_file', '-v', 'variants_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='sero_res', id='id', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))

    def test_arguments_short_options_surface_typing(self):
        actual = get_arguments().parse_args(['surface_typer', '-i', 'id', '-x', 'file1', '-y', 'file2', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='surface_typer', id='id', surface_inc='file1', surface_variants='file2', output='output_prefix'))

    def test_arguments_short_options_pbp_typing(self):
        actual = get_arguments().parse_args(['pbp_typer', '-i', 'id', '-p', 'pbp_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='pbp_typer', id='id', pbp_allele='pbp_file', output='output_prefix'))

    @patch('bin.combine_results.get_arguments')
    @patch('bin.combine_results.get_sero_res_contents')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.combine_results.get_content_with_id')
    def test_main_for_sero_res(self, mock_get_content_with_id, mock_write_output, mock_get_sero_res_contents, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "sero_res"
        mock_get_sero_res_contents.return_value = 'foobar1'
        mock_get_content_with_id.return_value = ANY
        main()
        self.assertEqual(mock_get_sero_res_contents.call_args_list, [call(args.id, args.sero, args.inc)])
        mock_write_output.assert_has_calls([
            call(ANY, args.output + "_sero_res_incidence.txt"),
            call(ANY, args.output + "_id_alleles_variants.txt"),
            call(ANY, args.output + "_id_variants.txt")
            ], any_order=False)
        mock_get_content_with_id.assert_has_calls([
            call(args.id, args.alleles),
            call(args.id, args.variants)
        ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.combine_results.get_content_with_id')
    def test_main_for_surface_typer(self, mock_get_content_with_id, mock_write_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "surface_typer"
        mock_get_content_with_id.return_value = ANY
        main()
        mock_write_output.assert_has_calls([
            call(ANY, args.output + "_surface_protein_variants.txt"),
            call(ANY, args.output + "_surface_protein_incidence.txt")
        ], any_order=False)
        mock_get_content_with_id.assert_has_calls([
            call(args.id, args.surface_inc),
            call(args.id, args.surface_variants)
        ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.combine_results.get_content_with_id')
    def test_main_for_pbp_typer(self, mock_get_content_with_id, mock_write_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "pbp_typer"
        mock_get_content_with_id.return_value = ANY
        main()
        mock_write_output.assert_has_calls([
            call(ANY, args.output + "_existing_PBP_allele.txt")
        ], any_order=False)
        mock_get_content_with_id.assert_has_calls([
            call(args.id, args.pbp_allele)
        ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.combine_results.get_content_with_id')
    def test_main_no_option(self, mock_get_content_with_id, mock_write_output, mock_get_arguments):
        main()
        mock_get_arguments.return_value.print_help.assert_called_once()
