from unittest import TestCase

import os
from bin.get_relative_abundance import *

class GetRelativeAbundance(TestCase):
    TEST_LANE_ID = 'test_lane1'
    TEST_LANE_ID2 = 'test_lane2'
    TEST_LANE_ID3 = 'test_lane3'
    TEST_DATA_DIR = 'tests/test_data'
    TEST_KRAKEN_REPORT = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/kraken.report'
    TEST_KRAKEN_REPORT_MULTIPLE_LEADING_WHITESPACE  = f'{TEST_DATA_DIR}/{TEST_LANE_ID2}/kraken.report'
    TEST_KRAKEN_REPORT_EMPTY = f'{TEST_DATA_DIR}/{TEST_LANE_ID3}/kraken.report'
    TEST_OUTPUT_FILE = f'tests/test_data/relative_abundance.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/relative_abundance2.tab'
    TEST_DEST_FILE = f'{TEST_DATA_DIR}/test_file_dest.txt'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'


    def test_search_rel_abund(self):
        actual = search_rel_abund(self.TEST_KRAKEN_REPORT, "Streptococcus agalactiae")

        self.assertEqual(92.38, actual)

    def test_get_relative_abundance(self):
        actual = get_relative_abundance(['tests/test_data/test_lane1', 'tests/test_data/test_lane2', 'tests/test_data/test_lane3'], [self.TEST_LANE_ID, self.TEST_LANE_ID2, self.TEST_LANE_ID3, "test_lane4"], "Streptococcus agalactiae")

        self.assertEqual([("test_lane1", 92.38), ("test_lane2", 2.38), ("test_lane3", None), ("test_lane4", None)], actual)

    def test_write_qc_status(self):
        write_qc_status([("lane1", 92.38), ("lane2", 2.38), ("lane3", None), ("lane4", None)], 70, [ "rel_abundance", "rel_abundance_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\trel_abundance\trel_abundance_status\nlane1\t92.38\tPASS\nlane2\t2.38\tFAIL\nlane3\t\t\nlane4\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--file_dest', 'kraken_reports', '--species', 'Streptococcus agalactiae', '--threshold', '70',
            '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='kraken_reports', species='Streptococcus agalactiae',
                         threshold=70.0, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-d', 'kraken_reports', '-s', 'Streptococcus agalactiae', '-t', '70',
            '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='kraken_reports', species='Streptococcus agalactiae',
                         threshold=70.0, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--file_dest', self.TEST_DEST_FILE, '--species', 'Streptococcus agalactiae', '--threshold', '70',
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\trel_abundance\trel_abundance_status\ntest_lane1\t92.38\tPASS\ntest_lane2\t2.38\tFAIL\ntest_lane3\t\t\ntest_lane4\t\t\n""")
