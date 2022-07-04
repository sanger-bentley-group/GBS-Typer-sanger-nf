import unittest
import argparse
from unittest.mock import patch, call, ANY
import pandas as pd
import numpy as np
import os

from bin.combine_results import read_header_json, get_content, create_model_df, merge_dfs, create_df, rename_columns, get_arguments, main
from lib.file_utils import FileUtils

class TestCombineResults(unittest.TestCase):

    TEST_LANE = '25292_2#85'
    TEST_HEADERS_FILE = 'headers.json'
    TEST_DATA_SEROTYPE = 'test_data/input/' + TEST_LANE + '_SeroType_Results.txt'
    TEST_DATA_RES_INCIDENCE = 'test_data/input/' + TEST_LANE + '_res_incidence.txt'
    TEST_DATA_RES_ALLELES = 'test_data/input/' + TEST_LANE + '_res_alleles.txt'
    TEST_DATA_RES_VARIANTS = 'test_data/input/' + TEST_LANE + '_res_gbs_variants.txt'
    TEST_DATA_PBP_ALLELE = 'test_data/input/' + TEST_LANE + '_GBS1A-1_PBP_existing_allele.txt'
    TEST_DATA_SURFACE_TYPER = 'test_data/input/' + TEST_LANE + '_surface_protein_incidence_sample.txt'
    TEST_DATA_EMPTY_SURFACE_TYPER = 'test_data/input/empty_surface_protein_incidence_sample.txt'
    TEST_DATA_MLST_ALLELIC_FREQUENCY = 'test_data/input/' + TEST_LANE + '__mlst__Streptococcus_agalactiae_MLST_alleles__results.txt'
    TEST_OUTPUT = "test_data/output/" + TEST_LANE + "_output.txt"

    header_dict = read_header_json(TEST_HEADERS_FILE)
    id_df = pd.DataFrame(TEST_LANE, columns=header_dict["id"], index = [0])

    def test_create_model_df_for_main_report(self):
        actual = create_model_df(self.header_dict["combine_all"].keys(), self.id_df)
        self.maxDiff = None
        expected = {
            'lane_id': {0: '25292_2#85'},
            'Serotype': {0: np.nan},
            'ST': {0: np.nan},
            'adhP': {0: np.nan},
            'pheS': {0: np.nan},
            'atr': {0: np.nan},
            'glnA': {0: np.nan},
            'sdhA': {0: np.nan},
            'glcK': {0: np.nan},
            'tkt': {0: np.nan},
            "AAC6APH2": {0: np.nan},
            'ANT6IA': {0: np.nan},
            "ANT6IA3KF864551": {0: np.nan},
            "APH3III": {0: np.nan},
            'CATPC194': {0: np.nan},
            'CATQ': {0: np.nan},
            'ERMA': {0: np.nan},
            'ERMB':	{0: np.nan},
            'ERMT':	{0: np.nan},
            'LNUB':	{0: np.nan},
            'LNUC': {0: np.nan},
            'LSAC':	{0: np.nan},
            'LSAE':	{0: np.nan},
            'MEFA':	{0: np.nan},
            'MSRD':	{0: np.nan},
            'TETB': {0: np.nan},
            'TETL': {0: np.nan},
            'TETM': {0: np.nan},
            'TETO': {0: np.nan},
            "TETO32O": {0: np.nan},
            "TETOW": {0: np.nan},
            "TETOW32O": {0: np.nan},
            "TETOW32OWO": {0: np.nan},
            "TETOWO": {0: np.nan},
            "TETS": {0: np.nan},
            "TETSM": {0: np.nan},
            "TETW32O": {0: np.nan},
            "TETW4FN396364": {0: np.nan},
            'ALP1': {0: np.nan},
            'ALP23': {0: np.nan},
            'ALPHA': {0: np.nan},
            'HVGA': {0: np.nan},
            'PI1': {0: np.nan},
            'PI2A1': {0: np.nan},
            'PI2A2': {0: np.nan},
            'PI2B': {0: np.nan},
            'RIB': {0: np.nan},
            'SRR1': {0: np.nan},
            'SRR2': {0: np.nan},
            "23S1_SNP": {0: np.nan},
            "23S3_SNP": {0: np.nan},
            'GYRA_SNP': {0: np.nan},
            'PARC_SNP': {0: np.nan}
        }
        self.assertEqual(actual.to_dict(), expected)

    def test_create_df_with_id_alleles(self):
        df_res_alleles = create_df(self.header_dict["res_alleles"], self.id_df, [self.TEST_DATA_RES_ALLELES])
        self.assertEqual(df_res_alleles.to_dict(), {
            "lane_id": {0: '25292_2#85'},
            "AG": {0: "aac(6')-30-aac(6')-Ib'[aac(6')-30-aac(6')-Ib'_1]"},
            'EC': {0: 'ERMB(ERMB-1):MEF(MEF-1):23S1:23S3'},
            "FQ": {0: 'PARC:GYRA-V1A,M2Q,G3W,K4W'},
            "OTHER": {0: "MsrD_MLS(MsrD_296):Ant6-Ia_AGly(Ant6-Ia_1633):Sat4A_AGly(Sat4A_586):Aph3-III_AGly(Aph3-III_268):msr(D)(msr(D)_2):ant(6)-Ia(ant(6)-Ia):aph(3')-III(aph(3')-III)"},
            'TET': {0: 'TETO(TETO-1)'}
        })

    def test_create_df_with_id_variants(self):
        df_gbs_res_variants = create_df(self.header_dict["gbs_res_variants"], self.id_df, [self.TEST_DATA_RES_VARIANTS])
        self.assertEqual(df_gbs_res_variants.to_dict(), {
            "lane_id": {0: '25292_2#85'},
            "23S1_SNP": {0: '*'},
            "23S3_SNP": {0: '*'},
            "GYRA_SNP": {0: 'V1A,M2Q,G3W,K4W'},
            "PARC_SNP": {0: '*'}
        })

    def test_create_df_with_sero_res(self):
        df_sero_res = create_df(self.header_dict["sero_res"], self.id_df, [self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE])
        self.maxDiff = None
        self.assertEqual(df_sero_res.to_dict(), {
            "lane_id": {0: '25292_2#85'},
            "Serotype": {0: 'III'},
            "AAC6APH2": {0: 'pos'},
            "ANT6IA": {0: 'neg'},
            'ANT6IA3KF864551': {0: 'neg'},
            "APH3III": {0: 'pos'},
            'CATPC194': {0: 'pos'},
            "CATQ": {0: 'neg'},
            "ERMA": {0: 'neg'},
            "ERMB": {0: 'neg'},
            "ERMT": {0: 'neg'},
            "LNUB": {0: 'neg'},
            "LNUC": {0: 'neg'},
            "LSAC": {0: 'neg'},
            "LSAE": {0: 'neg'},
            "MEFA": {0: 'neg'},
            "MSRD": {0: 'neg'},
            "TETB": {0: 'neg'},
            "TETL": {0: 'neg'},
            "TETM": {0: 'pos'},
            "TETO": {0: 'neg'},
            "TETO32O": {0: 'neg'},
            "TETOW": {0: 'neg'},
            "TETOW32O": {0: 'neg'},
            "TETOW32OWO": {0: 'neg'},
            "TETOWO": {0: 'neg'},
            "TETS": {0: 'neg'},
            "TETSM": {0: 'neg'},
            "TETW32O": {0: 'neg'},
            "TETW4FN396364": {0: 'neg'}
            })

    def test_create_df_for_all_content(self):
        df_combine_all = create_df(list(self.header_dict["combine_all"].keys()), self.id_df, [self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE, self.TEST_DATA_RES_VARIANTS, self.TEST_DATA_MLST_ALLELIC_FREQUENCY, self.TEST_DATA_SURFACE_TYPER])
        df_combine_all = df_combine_all.replace(to_replace=['+', '-'], value=['pos', 'neg'])
        df_combine_all = rename_columns(df_combine_all, self.header_dict["combine_all"], self.id_df)
        self.maxDiff = None
        self.assertEqual(list(df_combine_all.to_dict().keys()), ['lane_id', 'cps_type', 'ST', 'adhP', 'pheS', 'atr', 'glnA', 'sdhA', 'glcK', 'tkt', "aac(6')-aph(2'')", 'ant(6-Ia)', "aph(3'-III)", 'aadE', 'cat(pc194)', 'catQ', 'ermA', 'ermB', 'ermT', 'lnuB', 'lnuC', 'lsaC', 'lsaE', 'mefA', 'msrD', 'tetB', 'tetL', 'tetM', 'tetW', 'tetO', 'tetS', 'tetO32O', 'tetOW', 'tetOW32O', 'tetOW32OWO', 'tetOWO', 'tetSM', 'tetW32O', 'alp1', 'alp2/3', 'alpha', 'hvgA', 'PI1', 'PI2A1', 'PI2A2', 'PI2B', 'rib', 'srr1', 'srr2', '23S1_SNP', '23S3_SNP', 'gyrA_SNP', 'parC_SNP'])

    def test_create_df_for_empty_file(self):
        actual = create_df(self.header_dict["surface_inc"], self.id_df, [self.TEST_DATA_EMPTY_SURFACE_TYPER])
        self.assertEqual(actual.to_dict(),
        {
            'lane_id': {0: '25292_2#85'},
            'ALP1': {0: np.nan},
            'ALP23': {0: np.nan},
            'ALPHA': {0: np.nan},
            'HVGA': {0: np.nan},
            'PI1': {0: np.nan},
            'PI2A1': {0: np.nan},
            'PI2A2': {0: np.nan},
            'PI2B': {0: np.nan},
            'RIB': {0: np.nan},
            'SRR1': {0: np.nan},
            'SRR2': {0: np.nan}
        })

    def test_should_get_pbp_contents(self):
        actual = create_df(self.header_dict["pbp_allele"], self.id_df, [self.TEST_DATA_PBP_ALLELE])
        self.assertEqual(actual.to_dict(), {
            "lane_id": {0: "25292_2#85"},
            "Contig": {0: ".25292_2_85.9:39495-40455(+)"},
            "PBP_allele": {0: "2||GBS_1A"}
        })

    def test_write_pandas_output_for_all_content(self):
        df_combine_all = create_df(list(self.header_dict["combine_all"].keys()), self.id_df, [self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE, self.TEST_DATA_RES_VARIANTS, self.TEST_DATA_MLST_ALLELIC_FREQUENCY, self.TEST_DATA_SURFACE_TYPER])
        df_combine_all = df_combine_all.replace(to_replace=['+', '-'], value=['pos', 'neg'])
        df_combine_all = rename_columns(df_combine_all, self.header_dict["combine_all"], self.id_df)
        FileUtils.write_pandas_output(df_combine_all, self.TEST_OUTPUT)
        actual = pd.read_csv(self.TEST_OUTPUT, sep="\t")
        self.maxDiff = None
        self.assertEqual(list(actual.to_dict().keys()), ['lane_id', 'cps_type', 'ST', 'adhP', 'pheS', 'atr', 'glnA', 'sdhA', 'glcK', 'tkt', "aac(6')-aph(2'')", 'ant(6-Ia)', "aph(3'-III)", 'aadE', 'cat(pc194)', 'catQ', 'ermA', 'ermB', 'ermT', 'lnuB', 'lnuC', 'lsaC', 'lsaE', 'mefA', 'msrD', 'tetB', 'tetL', 'tetM', 'tetW', 'tetO', 'tetS', 'tetO32O', 'tetOW', 'tetOW32O', 'tetOW32OWO', 'tetOWO', 'tetSM', 'tetW32O', 'alp1', 'alp2/3', 'alpha', 'hvgA', 'PI1', 'PI2A1', 'PI2A2', 'PI2B', 'rib', 'srr1', 'srr2', '23S1_SNP', '23S3_SNP', 'gyrA_SNP', 'parC_SNP'])
        os.remove(self.TEST_OUTPUT)

    def test_write_pandas_output_for_all_content_and_empty_surface_protein_file(self):
        df_combine_all = create_df(list(self.header_dict["combine_all"].keys()), self.id_df, [self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE, self.TEST_DATA_RES_VARIANTS, self.TEST_DATA_MLST_ALLELIC_FREQUENCY, self.TEST_DATA_EMPTY_SURFACE_TYPER])
        df_combine_all = df_combine_all.replace(to_replace=['+', '-'], value=['pos', 'neg'])
        df_combine_all = rename_columns(df_combine_all, self.header_dict["combine_all"], self.id_df)
        FileUtils.write_pandas_output(df_combine_all, self.TEST_OUTPUT)
        actual = pd.read_csv(self.TEST_OUTPUT, sep="\t")
        self.maxDiff = None
        self.assertEqual(list(actual.to_dict().keys()), ['lane_id', 'cps_type', 'ST', 'adhP', 'pheS', 'atr', 'glnA', 'sdhA', 'glcK', 'tkt', "aac(6')-aph(2'')", 'ant(6-Ia)', "aph(3'-III)", 'aadE', 'cat(pc194)', 'catQ', 'ermA', 'ermB', 'ermT', 'lnuB', 'lnuC', 'lsaC', 'lsaE', 'mefA', 'msrD', 'tetB', 'tetL', 'tetM', 'tetW', 'tetO', 'tetS','tetO32O', 'tetOW', 'tetOW32O', 'tetOW32OWO','tetOWO','tetSM','tetW32O','alp1', 'alp2/3', 'alpha', 'hvgA', 'PI1', 'PI2A1', 'PI2A2', 'PI2B', 'rib', 'srr1', 'srr2', '23S1_SNP', '23S3_SNP', 'gyrA_SNP', 'parC_SNP'])
        os.remove(self.TEST_OUTPUT)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_pandas_output')
    @patch('bin.combine_results.create_df')
    @patch('bin.combine_results.read_header_json')
    def test_main_for_sero_res(self, mock_read_header_json, mock_create_df, mock_write_pandas_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "sero_res"
        header_dict = self.header_dict
        mock_read_header_json.return_value = header_dict
        mock_create_df.return_value = 'foobar1'

        main()

        mock_create_df.assert_has_calls([
            call(header_dict['sero_res'], ANY, [args.sero, args.inc]),
            call(header_dict['res_alleles'], ANY, [args.alleles]),
            call(header_dict['gbs_res_variants'], ANY, [args.variants])
            ], any_order=False)
        mock_write_pandas_output.assert_has_calls([
            call('foobar1', args.output + "_sero_res_incidence.txt"),
            call('foobar1', args.output + "_id_alleles_variants.txt"),
            call('foobar1', args.output + "_id_variants.txt")
            ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_pandas_output')
    @patch('bin.combine_results.create_df')
    @patch('bin.combine_results.read_header_json')
    def test_main_for_surface_typer(self, mock_read_header_json, mock_create_df, mock_write_pandas_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "surface_typer"
        header_dict = self.header_dict
        mock_read_header_json.return_value = header_dict
        mock_create_df.return_value = 'foobar1'

        main()

        mock_create_df.assert_has_calls([
            call(header_dict['surface_inc'], ANY, [args.surface_inc]),
            call(header_dict['surface_variants'], ANY, [args.surface_variants])
            ], any_order=False)
        mock_write_pandas_output.assert_has_calls([
            call('foobar1', args.output + "_surface_protein_incidence.txt"),
            call('foobar1', args.output + "_surface_protein_variants.txt")
            ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_pandas_output')
    @patch('bin.combine_results.rename_columns')
    @patch('bin.combine_results.create_df')
    @patch('bin.combine_results.read_header_json')
    def test_main_for_combine_all(self, mock_read_header_json, mock_create_df, mock_rename_columns, mock_write_pandas_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "combine_all"
        header_dict = self.header_dict
        mock_read_header_json.return_value = header_dict
        mock_create_df.return_value = 'foobar1'
        mock_rename_columns.return_value = 'renamed'

        main()

        mock_create_df.assert_has_calls([
            call(list(header_dict["combine_all"].keys()), ANY, [args.sero, args.inc, args.variants, args.mlst, args.surface_inc])], any_order=False)
        mock_rename_columns.assert_has_calls([
            call('foobar1', header_dict['combine_all'], ANY)
        ])
        mock_write_pandas_output.assert_has_calls([
            call('renamed', args.output + "_id_combined_output.txt")
        ], any_order=False)

    @patch('bin.combine_results.get_arguments')
    @patch('lib.file_utils.FileUtils.write_pandas_output')
    @patch('bin.combine_results.create_df')
    @patch('bin.combine_results.read_header_json')
    def test_main_for_pbp_typer(self, mock_read_header_json, mock_create_df, mock_write_pandas_output, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.which = "pbp_typer"
        header_dict = self.header_dict
        mock_read_header_json.return_value = header_dict
        mock_create_df.return_value = 'foobar1'

        main()

        mock_create_df.assert_has_calls([
            call(header_dict['pbp_allele'], ANY, [args.pbp_allele])], any_order=False)
        mock_write_pandas_output.assert_has_calls([
            call('foobar1', args.output + "_existing_PBP_allele.txt")
            ], any_order=False)

    def test_arguments_sero_res(self):
        actual = get_arguments().parse_args(['sero_res', '--id', 'id', '--headers', 'header_file', '--serotyper_results', 'sero_file', '--res_incidence_results', 'res_file',
                                            '--res_alleles_results', 'allele_file', '--res_variants_results', 'variants_file', '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='sero_res', id='id', headers='header_file', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))

    def test_arguments_surface_typing(self):
        actual = get_arguments().parse_args(['surface_typer', '--id', 'id', '--headers', 'header_file', '--surface_incidence_results', 'file1', '--surface_variants_results', 'file2',
                                             '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='surface_typer', id='id', headers='header_file', surface_inc='file1', surface_variants='file2', output='output_prefix'))

    def test_arguments_pbp_typing(self):
        actual = get_arguments().parse_args(['pbp_typer', '--id', 'id', '--headers', 'header_file', '--pbp_existing_allele_results', 'pbp_file', '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='pbp_typer', id='id', headers='header_file', pbp_allele='pbp_file', output='output_prefix'))

    def test_arguments_short_options_sero_res(self):
        actual = get_arguments().parse_args(['sero_res', '-i', 'id', '-t', 'header_file', '-s', 'sero_file', '-r', 'res_file', '-a', 'allele_file', '-v', 'variants_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='sero_res', id='id', headers='header_file', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))

    def test_arguments_short_options_surface_typing(self):
        actual = get_arguments().parse_args(['surface_typer', '-i', 'id', '-t', 'header_file', '-x', 'file1', '-y', 'file2', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='surface_typer', id='id', headers='header_file', surface_inc='file1', surface_variants='file2', output='output_prefix'))

    def test_arguments_short_options_pbp_typing(self):
        actual = get_arguments().parse_args(['pbp_typer', '-i', 'id', '-t', 'header_file', '-p', 'pbp_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='pbp_typer', id='id', headers='header_file', pbp_allele='pbp_file', output='output_prefix'))

    def test_arguments_short_options_combine_all(self):
        actual = get_arguments().parse_args(['combine_all', '-i', 'id', '-t', 'header_file', '-s', 'sero_file', '-r', 'res_file', '-v', 'variants_file', '-m', 'mlst_file', '-x', 'surface_typer_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(which='combine_all', id='id', headers='header_file', sero='sero_file', inc='res_file', variants='variants_file', mlst='mlst_file', surface_inc='surface_typer_file', output='output_prefix'))
