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
