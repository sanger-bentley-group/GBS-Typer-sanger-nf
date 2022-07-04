from unittest import TestCase

import os
from bin.get_genome_length import *

class GetGenomeLength(TestCase):
    TEST_LANE_ID = 'test_lane1'
    TEST_DATA_DIR = 'tests/test_data'
    TEST_VELVET_FASTA = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/velvet_assembly/contigs.fa'
    TEST_OUTPUT_FILE = f'tests/test_data/genome_length.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/genome_length2.tab'
    TEST_DEST_FILE = f'{TEST_DATA_DIR}/test_file_dest.txt'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'


    def test_count_nucl(self):
        actual = count_nucl(self.TEST_VELVET_FASTA)

        self.assertEqual(2000, actual)

    def test_get_genome_length(self):
        actual = get_genome_length(['tests/test_data/test_lane1', 'tests/test_data/test_lane2', 'tests/test_data/test_lane3'], [self.TEST_LANE_ID, "test_lane2", "test_lane3", "test_lane4"], "velvet")

        self.assertEqual([("test_lane1", 2000), ("test_lane2", None), ("test_lane3", None), ("test_lane4", None)], actual)

    def test_write_qc_status(self):
        write_qc_status([("test_lane1", 2000), ("test_lane2", None), ("test_lane3", None)], 1999, 2001, [ "genome_len", "genome_len_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\tgenome_len\tgenome_len_status\ntest_lane1\t2000\tPASS\ntest_lane2\t\t\ntest_lane3\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--file_dest', 'file_dest', '--assembler', 'spades', '--lower_threshold', '290000',
            '--higher_threshold', '1450000', '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         lower_threshold=290000, higher_threshold=1450000, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-d', 'file_dest', '-a', 'spades', '-lt', '290000',
            '-ht', '1450000', '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         lower_threshold=290000, higher_threshold=1450000, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--file_dest', self.TEST_DEST_FILE, '--assembler', 'velvet', '--lower_threshold', '290000', '--higher_threshold', '1450000',
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\tgenome_len\tgenome_len_status\ntest_lane1\t2000\tFAIL\ntest_lane2\t\t\ntest_lane3\t\t\ntest_lane4\t\t\n""")
