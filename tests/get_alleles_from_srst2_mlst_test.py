import argparse
import os
import unittest

from bin.get_alleles_from_srst2_mlst import get_mismatch_and_depth, get_new_and_existing_alleles, write_alleles_file, get_arguments
class TestProcessResults(unittest.TestCase):
    TEST_OUTPUT_PREFIX = 'test'
    TEST_FILE = 'tests/test_data/input/test__mlst__Streptococcus_agalactiae_MLST_alleles__results.txt'
    TEST_OUT1 = 'tests/test_data/output/test_mlst_alleles.txt'
    TEST_OUT2 = 'tests/test_data/output/test_mlst_alleles2.txt'

    def test_get_mismatch_and_depth(self):
        actual = get_mismatch_and_depth(self.TEST_FILE)
        self.assertEqual(actual, ('adhP_1/1snp', 173.614142857, 'ST-1'))

    def test_get_new_and_existing_alleles(self):
        actual = get_new_and_existing_alleles(('adhP_1/1snp', 173.614142857, 'ST-1'), 30, self.TEST_OUTPUT_PREFIX)
        f = open(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt', "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """Alleles found\nadhP_1\n""")
        os.remove(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt')

    def test_get_new_and_existing_alleles_low_depth(self):
        actual = get_new_and_existing_alleles(('adhP_1/1snp', 29.99, 'ST-1'), 30, self.TEST_OUTPUT_PREFIX)
        f = open(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt', "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """test: No new MLST alleles were found with sufficient read depth above 30.\n""")
        os.remove(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt')

    def test_get_new_and_existing_alleles_multi_alleles(self):
        actual = get_new_and_existing_alleles(('adhP_1/1snp;pheS_1/1snp', 173.614142857, 'ST-1'), 30, self.TEST_OUTPUT_PREFIX)
        f = open(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt', "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """Alleles found\nadhP_1\npheS_1\n""")
        os.remove(self.TEST_OUTPUT_PREFIX + '_new_mlst_alleles.txt')

    def test_get_new_and_existing_alleles_no_mismatches(self):
        actual = get_new_and_existing_alleles(('0', 173.614142857, 'ST-1'), 30, self.TEST_OUTPUT_PREFIX)
        f = open(self.TEST_OUTPUT_PREFIX + '_existing_mlst_alleles.txt', "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """ID\tST\ntest\tST-1\n""")
        os.remove(self.TEST_OUTPUT_PREFIX + '_existing_mlst_alleles.txt')

    def test_alleles_file(self):
        write_alleles_file(['Alleles found', 'adhP_1', 'pheS_1'], self.TEST_OUT1)
        f = open(self.TEST_OUT1, "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """Alleles found\nadhP_1\npheS_1\n""")

    def test_alleles_file_without_alleles(self):
        write_alleles_file(['No new MLST alleles were found.'], self.TEST_OUT2)
        f = open(self.TEST_OUT2, "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """No new MLST alleles were found.\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--mlst_results_file', 'mlst_file', '--min_read_depth', '30',
            '--output_prefix', 'out'])
        self.assertEqual(actual,
                         argparse.Namespace(mlst='mlst_file', min_depth=30, output='out'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-m', 'mlst_file', '-d', '30', '-o', 'out'])
        self.assertEqual(actual,
                         argparse.Namespace(mlst='mlst_file', min_depth=30, output='out'))
