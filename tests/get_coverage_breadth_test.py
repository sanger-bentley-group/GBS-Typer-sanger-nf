from unittest import TestCase

import os
from bin.get_coverage_breadth import *

class GetCoverageBreadth(TestCase):
    TEST_DATA_DIR = 'tests/test_data'
    TEST_PF_STATS = f'{TEST_DATA_DIR}/test.txt.pathfind_stats.csv'
    TEST_OUTPUT_FILE = f'tests/test_data/breadth_of_coverage.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/breadth_of_coverage2.tab'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'

    lane_ids_cov_breadth = defaultdict(lambda: None)
    lane_ids_cov_breadth['test_lane1'] = 66.83
    lane_ids_cov_breadth['test_lane2'] = 86.46
    lane_ids_cov_breadth['test_lane3'] = 86.7


    def test_get_coverage_breadth(self):
        actual = get_coverage_breadth(self.TEST_PF_STATS)

        self.assertEqual(self.lane_ids_cov_breadth, actual)

    def test_write_qc_status(self):
        write_qc_status(['test_lane1', 'test_lane2', 'test_lane3', 'test_lane4'], self.lane_ids_cov_breadth, 70, [ "cov_breadth", "cov_breadth_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\tcov_breadth\tcov_breadth_status\ntest_lane1\t66.83\tFAIL\ntest_lane2\t86.46\tPASS\ntest_lane3\t86.7\tPASS\ntest_lane4\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--qc_stats', 'qc_stats', '--threshold', '70',
            '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', qc_stats='qc_stats',
                         threshold=70, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-q', 'qc_stats', '-t', '70',
            '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', qc_stats='qc_stats',
                         threshold=70, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--qc_stats', self.TEST_PF_STATS, '--threshold', '70',
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\tcov_breadth\tcov_breadth_status\ntest_lane1\t66.83\tFAIL\ntest_lane2\t86.46\tPASS\ntest_lane3\t86.7\tPASS\ntest_lane4\t\t\n""")
