from unittest import TestCase

import os
from bin.get_gc_content import *

class GetGCContent(TestCase):
    TEST_LANE_ID = 'test_lane1'
    TEST_DATA_DIR = 'tests/test_data'
    TEST_SPADES_FASTA = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/spades_assembly/contigs.fa'
    TEST_VELVET_FASTA = f'{TEST_DATA_DIR}/{TEST_LANE_ID}/velvet_assembly/contigs.fa'
    TEST_OUTPUT_FILE = f'tests/test_data/gc_content.tab'
    TEST_OUTPUT_FILE2 = f'tests/test_data/gc_content2.tab'
    TEST_DEST_FILE = f'{TEST_DATA_DIR}/test_file_dest.txt'
    TEST_LANE_IDS = f'{TEST_DATA_DIR}/test_lane_ids.txt'
    TEST_HEADERS_FILE = f'{TEST_DATA_DIR}/test_headers.json'


    def test_gc_content_spades(self):
        actual = gc_content(self.TEST_SPADES_FASTA)

        self.assertEqual(50.0, actual)

    def test_gc_content_velvet(self):
        actual = gc_content(self.TEST_VELVET_FASTA)

        self.assertEqual(50.0, actual)

    def test_gc_content(self):
        actual = get_gc_content(['tests/test_data/test_lane1', 'tests/test_data/test_lane2', 'tests/test_data/test_lane3'], [self.TEST_LANE_ID, "test_lane2", "test_lane3", "test_lane4"], "spades")

        self.assertEqual([("test_lane1", 50.0), ("test_lane2", None), ("test_lane3", None), ("test_lane4", None)], actual)

    def test_write_qc_status(self):
        write_qc_status([("test_lane1", 50.0), ("test_lane2", 33.3), ("test_lane3", None)], 32, 38, [ "gc_content", "gc_content_status"], self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual, """lane_id\tgc_content\tgc_content_status\ntest_lane1\t50.0\tFAIL\ntest_lane2\t33.3\tPASS\ntest_lane3\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--lane_ids', 'lane_file', '--file_dest', 'file_dest', '--assembler', 'spades', '--lower_threshold', '32',
            '--higher_threshold', '38', '--headers', 'header_json', '--output_file', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         lower_threshold=32, higher_threshold=38, headers='header_json', output_file='output_tab_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-l', 'lane_file', '-d', 'file_dest', '-a', 'spades', '-lt', '32',
            '-ht', '38', '-j', 'header_json', '-o', 'output_tab_file'])
        self.assertEqual(actual,
                         argparse.Namespace(lane_ids='lane_file', file_dest='file_dest', assembler='spades',
                         lower_threshold=32, higher_threshold=38, headers='header_json', output_file='output_tab_file'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--lane_ids', self.TEST_LANE_IDS, '--file_dest', self.TEST_DEST_FILE, '--assembler', 'spades', '--lower_threshold', '32', '--higher_threshold', '38',
            '--headers', self.TEST_HEADERS_FILE, '--output_file', self.TEST_OUTPUT_FILE2])

        main(args)

        file = open(self.TEST_OUTPUT_FILE2, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE2)

        self.assertEqual(actual, """lane_id\tgc_content\tgc_content_status\ntest_lane1\t50.0\tFAIL\ntest_lane2\t\t\ntest_lane3\t\t\ntest_lane4\t\t\n""")
