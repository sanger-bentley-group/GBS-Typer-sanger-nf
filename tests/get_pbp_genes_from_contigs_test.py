import argparse
import unittest
from unittest.mock import patch, call, ANY

from bin.get_pbp_genes_from_contigs import BlastData, SeqData, FragmentData, get_arguments, check_arguments, main

class TestGetPBPGenesFromContigs(unittest.TestCase):
    TEST_BLAST_DATA = 'test_data/input/test_blast_blactam.out'
    TEST_SEQ_DATA = 'test_data/input/test_GBS_bLactam_Ref.fasta'
    TEST_OUTPUT_PREFIX = 'test_data/output/TEST_'


    def test_read_blast_out(self):
        """
        Test output of blast
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        params_list = [
            ('GBS1A-1', 0, ['.26077_6_118.11', '100.000', '960', '0', '0', '1', '960', '39459', '40418', '0.0', '1773']),
            ('GBS2B-1', 1, ['.26077_6_118.2', '100.000', '16', '0', '0', '534', '549', '101982', '101967', '1.3', '30.7']),
            ('GBS2X-1', 2, ['.26077_6_118.10', '100.000', '15', '0', '0', '611', '625', '8263', '8277', '4.5', '28.8'])
        ]
        for param in params_list:
            actual = blast_data_to_process.get_data()[param[0]][param[1]]
            self.assertEqual(actual, param[2])


    def test_get_best_blast_hit(self):
        """
        Test the best blast hit of each beta lactam
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        params_list = [
            ('GBS1A-1', ['.26077_6_118.11', '100.000', '960', '0', '0', '1', '960', '39459', '40418', '0.0', '1773']),
            ('GBS2B-1', ['.26077_6_118.2', '100.000', '1065', '0', '0', '1', '1065', '185772', '186836', '0.0',	'1967']),
            ('GBS2X-1', ['.26077_6_118.11', '100.000', '1038', '0', '0', '1', '1038', '52297', '51260',	'0.0', '1917'])
        ]
        for param in params_list:
            self.assertEqual(blast_data_to_process.get_best_hit()[param[0]], param[1])


    def test_get_start_end_positions(self):
        """
        Test getting the start and end blactam positions in the contigs
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        seq_lengths = seq_data_to_analyse.calculate_seq_length()
        best_blast_hits = blast_data_to_process.get_best_hit()
        fragment_data = FragmentData()
        fragment_data.get_start_end_positions(best_blast_hits, seq_lengths, 0.5, 0.5)

        self.assertEqual(fragment_data.get_data(), {'GBS1A-1': ('.26077_6_118.11', '39458', '40418', 'forward', '1', '+'),
                                                    'GBS2B-1': ('.26077_6_118.2', '185771', '186836', 'forward', '1', '+'),
                                                    'GBS2X-1': ('.26077_6_118.11', '51259', '52297', 'reverse', '1', '-')})


    def test_write_start_end_positions(self):
        """
        Test output file contents containing start and end positions
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        seq_lengths = seq_data_to_analyse.calculate_seq_length()
        best_blast_hits = blast_data_to_process.get_best_hit()
        fragment_data = FragmentData()
        fragment_data.get_start_end_positions(best_blast_hits, seq_lengths, 0.5, 0.5)

        params_list = [
            ('GBS1A-1', ['.26077_6_118.11\t39458\t40418\tforward\t1\t+\n']),
            ('GBS2B-1', ['.26077_6_118.2\t185771\t186836\tforward\t1\t+\n']),
            ('GBS2X-1', ['.26077_6_118.11\t51259\t52297\treverse\t1\t-\n'])
        ]
        for param in params_list:
            fragment_data_location = 'test_data/output/' + param[0] + '.bed'
            fragment_data.write_start_end_positions('test_data/output/')
            fo = open(fragment_data_location, 'r')

            self.assertEqual(fo.readlines(), param[1])


    def test_read_seq_data(self):
        """
        Test output of sequence data
        """
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        actual = seq_data_to_analyse.get_data()

        self.assertEqual(actual['GBS1A-1'], 'GACATCTACAACAGTGACACTTACATCGCTTATCCAAACAATGAATTACAAATAGCATCTACCATCATGGATGCGACTAATGGTAAAGTCATTGCACAATTAGGCGGGCGTCATCAGAATGAAAATATTTCATTTGGGACAAATCAATCTGTCTTAACAGACCGCGATTGGGGTTCTACAATGAAACCTATCTCAGCTTATGCACCTGCTATTGATAGTGGTGTCTATAATTCAACAGGTCAATCATTAAACGACTCAGTTTACTACTGGCCTGGTACTTCTACTCAACTATATGACTGGGATCGTCAATATATGGGTTGGATGAGTATGCAGACCGCTATTCAACAATCACGTAACGTCCCTGCTGTCAGAGCACTTGAAGCCGCTGGATTAGACGAAGCAAAATCTTTCCTTGAAAAATTAGGCATATACTATCCAGAAATGAACTATTCAAATGCTATTTCAAGTAACAACAGTAGCAGTGATGCAAAATATGGTGCAAGTAGTGAGAAAATGGCAGCGGCTTACTCGGCTTTTGCAAACGGCGGAACTTACTATAAACCGCAATATGTTAATAAAATTGAATTTAGCGATGGAACCAATGATACTTATGCAGCGTCTGGTAGCCGTGCGATGAAAGAGACTACTGCCTACATGATGACGGATATGCTGAAAACAGTACTAACATTTGGTACTGGTACTAAAGCAGCTATCCCTGGTGTTGCACAAGCTGGTAAGACTGGTACTTCCAACTATACGGAAGATGAGTTAGCTAAAATTGAAGCAACTACTGGTATCTACAATAGCGCCGTTGGTACAATGGCTCCTGATGAAAACTTTGTCGGCTATACTTCTAAGTACACAATGGCAATTTGGACTGGTTATAAAAATCGCCTTACACCACTTTATGGTAGCCAACTGGATATTGCTACTGAGGTTTATCGTGCAATGATGTCCTAC')


    def test_calculate_seq_length(self):
        """
        Test the sequence lengths output
        """
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        actual = seq_data_to_analyse.calculate_seq_length()
        params_list = [
            ('GBS1A-1', 960),
            ('GBS2B-1', 1065),
            ('GBS2X-1', 1038)
        ]

        for param in params_list:
            self.assertEqual(actual[param[0]], param[1])


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--blast_out_file', 'blast_out_file', '--query_fasta', 'fasta_file',
            '--frac_align_len_threshold', '0.6', '--frac_identity_threshold', '0.6',
            '--output_prefix', 'out_prefix'])
        self.assertEqual(actual, argparse.Namespace(blast_out='blast_out_file',
        fasta_qu='fasta_file', frac_align=0.6, frac_ident=0.6, output='out_prefix'))


    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-b', 'blast_out_file', '-f', 'fasta_file',
            '-fa', '0.5', '-fi', '0.5', '-o', 'out_prefix'])
        self.assertEqual(actual, argparse.Namespace(blast_out='blast_out_file',
        fasta_qu='fasta_file', frac_align=0.5, frac_ident=0.5, output='out_prefix'))


    def test_check_arguments_frac_align(self):
        params_list = [
            (1.1, 0.5),
            (-0.1, 0.5)
        ]
        for param in params_list:
            args = argparse.Namespace(blast_out='blast_out_file',
            fasta_qu='fasta_file', frac_align=param[0], frac_ident=param[1], output='bed_file')

            with self.assertRaises(Exception) as exp:
                check_arguments(args)
            self.assertEqual(str(exp.exception), "Invalid frac_align_len_threshold value. Value must be between 0 and 1.")


    def test_check_arguments_frac_ident(self):
        params_list = [
            (0.5, 1.1),
            (0.5, -0.1)
        ]
        for param in params_list:
            args = argparse.Namespace(blast_out='blast_out_file',
            fasta_qu='fasta_file', frac_align=param[0], frac_ident=param[1], output='bed_file')

            with self.assertRaises(Exception) as exp:
                check_arguments(args)
            self.assertEqual(str(exp.exception), "Invalid frac_identity_threshold value. Value must be between 0 and 1.")


    @patch('bin.get_pbp_genes_from_contigs.get_arguments')
    def test_main(self, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.blast_out = self.TEST_BLAST_DATA
        args.fasta_qu = self.TEST_SEQ_DATA
        args.frac_align = 0.5
        args.frac_ident = 0.5
        args.output = self.TEST_OUTPUT_PREFIX
        main()
        fo = open('test_data/output/TEST_GBS1A-1.bed', 'r')

        self.assertEqual(fo.readlines(), ['.26077_6_118.11\t39458\t40418\tforward\t1\t+\n'])
