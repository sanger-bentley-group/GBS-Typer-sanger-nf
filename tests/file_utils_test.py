import unittest
import os
from lib.file_utils import FileUtils


class TestFileUtils(unittest.TestCase):
    """Unit test class for the file_utils module"""

    TEST_LANE = "26189_8#5"
    TEST_OUTPUT = "test_data/" + TEST_LANE + "_output.txt"

    def test_write_output(self):
        FileUtils.write_output('foobar', self.TEST_OUTPUT)
        f = open(self.TEST_OUTPUT, "r")
        actual = "".join(f.readlines())
        os.remove(self.TEST_OUTPUT)
        self.assertEqual(actual, """foobar""")

    def test_create_output_contents(self):
        final_dict = {'B_ITEM': 'pos', '1ITEM': 'neg', 'A_ITEM': 'neg'}
        actual = FileUtils.create_output_contents(final_dict)
        self.assertEqual(actual, '1ITEM\tA_ITEM\tB_ITEM\nneg\tneg\tpos\n')

    def test_read_delimited_id_file_with_hdrs(self):
        headers, contents = FileUtils.read_delimited_id_file_with_hdrs(
            '\t', 'test_data/input/delimited_id_file_with_hdrs.txt', 4)
        print(contents)
        self.assertEqual(3, len(headers))
        self.assertTrue('EC' in headers)
        self.assertTrue('FQ' in headers)
        self.assertTrue('OTHER' in headers)
        self.assertEqual(3, len(contents)) # num of lines
        self.assertEqual(3, len(contents['25292_2#85']))
        self.assertEqual(3, len(contents['25292_2#86']))
        self.assertEqual(3, len(contents['25292_2#87']))
        self.assertEqual('field1_1', contents['25292_2#85']['EC'])
        self.assertEqual('field1_2', contents['25292_2#85']['FQ'])
        self.assertEqual('field1_3', contents['25292_2#85']['OTHER'])
        self.assertEqual('field2_1', contents['25292_2#86']['EC'])
        self.assertEqual('field2_2', contents['25292_2#86']['FQ'])
        self.assertEqual('field2_3', contents['25292_2#86']['OTHER'])
        self.assertEqual('field3_1', contents['25292_2#87']['EC'])
        self.assertEqual('field3_2', contents['25292_2#87']['FQ'])
        self.assertEqual('field3_3', contents['25292_2#87']['OTHER'])

    def test_read_delimited_id_file_with_hdrs_bad_format(self):
        with self.assertRaises(RuntimeError):
            FileUtils.read_delimited_id_file_with_hdrs('\t', 'test_data/input/delimited_id_file_with_hdrs.txt', 10)
