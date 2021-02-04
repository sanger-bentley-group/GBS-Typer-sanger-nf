import argparse
import io
import unittest

from bin.process_serotyper_results import write_line, make_gene_dict, get_arguments



class TestProcessResults(unittest.TestCase):

    def test_should_make_gene_dict(self):
        actual = make_gene_dict('test_data/input/SERO_26237_7#5__fullgenes__GBS_seroT_Gene-DB_Final__results.txt', 10)
        self.assertEqual(actual, {'III': ['III', 'III-1', '100.0', '267.595', '1snp', '', '0.581', '172', '0.011', '4', '4']}
                         )

    def setUp(self) -> None:
        self.test_stream = io.StringIO()

    def test_should_write_line_imperfect(self):
        write_line('III', {'III': ['III', 'III-1', '100.0', '267.595', '1snp', '', '0.581', '172', '0.011', '4', '4']}, self.test_stream)
        actual = self.extract_written_content()
        self.assertEqual(actual, """III-1\tIII=imperfect\tIII\t267.595\n""")

    def test_should_write_line_identical(self):
        write_line('III', {'III': ['III', 'III-1', '100.0', '267.595', '', '', '0.581', '172', '0.011', '4', '4']}, self.test_stream)
        actual = self.extract_written_content()
        self.assertEqual(actual, """III-1\tIII=identical\tIII\t267.595\n""")

    def extract_written_content(self):
        self.test_stream.seek(0)
        return "".join(self.test_stream.readlines())

    def test_arguments(self):
        actual = get_arguments().parse_args(['--srst2_output', 'srst2_output_name', '--sero_db', 'sero_db', '--output', 'outfile', '--min_read_depth', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=30.0))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(['-s', 'srst2_output_name', '-b', 'sero_db', '-o', 'outfile', '-d', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=30.0))

    def test_arguments_without_depth_threshold(self):
        actual = get_arguments().parse_args(['-s', 'srst2_output_name', '-b', 'sero_db', '-o', 'outfile'])
        self.assertEqual(actual,
                        argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=30.0))
