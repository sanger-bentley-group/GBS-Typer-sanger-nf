from unittest import TestCase

import os
from bin.get_het_snps import *

class GetPercentageHETSNPs(TestCase):
    TEST_DATA_DIR = 'tests/test_data'
    TEST_STATS_DIR = f'{TEST_DATA_DIR}/het_stats'
    TEST_OUTPUT_FILE = f'tests/test_data/het_snps.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/het_snps2.tab'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'

    lane_ids_het_snps = defaultdict(lambda: None)
    lane_ids_het_snps['test_lane1'] = 4
    lane_ids_het_snps['test_lane2'] = 20
    lane_ids_het_snps['test_lane3'] = 21


    def test_get_het_snps(self):
        actual = get_het_snps(self.TEST_STATS_DIR)

        self.assertEqual(self.lane_ids_het_snps, actual)

    def test_write_qc_status(self):
        write_qc_status(['test_lane1', 'test_lane2', 'test_lane3', 'test_lane4'], self.lane_ids_het_snps, 20, [ "HET_SNPs", "HET_SNPs_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\tHET_SNPs\tHET_SNPs_status\ntest_lane1\t4\tPASS\ntest_lane2\t20\tPASS\ntest_lane3\t21\tFAIL\ntest_lane4\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--threshold', '20', '--data_dir', './',
            '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', data_dir='./',
                         threshold=20, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-t', '20', '-d', './',
            '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', data_dir='./',
                         threshold=20, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--threshold', '20', '--data_dir', self.TEST_STATS_DIR,
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\tHET_SNPs\tHET_SNPs_status\ntest_lane1\t4\tPASS\ntest_lane2\t20\tPASS\ntest_lane3\t21\tFAIL\ntest_lane4\t\t\n""")
