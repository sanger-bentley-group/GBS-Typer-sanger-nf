from unittest import TestCase

import os
from bin.get_no_contigs import *

class GetNoContigs(TestCase):
    TEST_LANE_ID = 'test_lane1'
    TEST_DATA_DIR = 'tests/test_data'
    TEST_SPADES_FASTA = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/spades_assembly/contigs.fa'
    TEST_VELVET_FASTA = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/velvet_assembly/contigs.fa'
    TEST_OUTPUT_FILE = f'tests/test_data/contig_number.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/contig_number2.tab'
    TEST_DEST_FILE = f'{TEST_DATA_DIR}/test_file_dest.txt'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'


    def test_count_contigs_spades(self):
        actual = count_contigs(self.TEST_SPADES_FASTA)

        self.assertEqual(1, actual)

    def test_count_contigs_velvet(self):
        actual = count_contigs(self.TEST_VELVET_FASTA)

        self.assertEqual(500, actual)

    def test_get_contig_number(self):
        actual = get_contig_number(['tests/test_data/test_lane1', 'tests/test_data/test_lane2', 'tests/test_data/test_lane3'], [self.TEST_LANE_ID, "test_lane2", "test_lane3", "test_lane4"], "spades")

        self.assertEqual([("test_lane1", 1), ("test_lane2", None), ("test_lane3", None), ("test_lane4", None)], actual)

    def test_write_qc_status(self):
        write_qc_status([("test_lane1", 1), ("test_lane2", 500), ("test_lane3", None)], 500, [ "contig_no", "contig_no_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\tcontig_no\tcontig_no_status\ntest_lane1\t1\tPASS\ntest_lane2\t500\tFAIL\ntest_lane3\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--file_dest', 'file_dest', '--assembler', 'spades', '--threshold', '500',
            '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         threshold=500, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-d', 'file_dest', '-a', 'spades', '-t', '500',
            '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         threshold=500, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--file_dest', self.TEST_DEST_FILE, '--assembler', 'spades', '--threshold', '500',
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\tcontig_no\tcontig_no_status\ntest_lane1\t1\tPASS\ntest_lane2\t\t\ntest_lane3\t\t\ntest_lane4\t\t\n""")
