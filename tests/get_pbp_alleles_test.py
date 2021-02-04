import argparse
import unittest
from unittest.mock import patch, call, ANY

from bin.get_pbp_alleles import BlastData, SeqData, get_identical_allele, get_imperfect_allele, write_content, get_arguments, main

class TestGetPBPAlleles(unittest.TestCase):
    TEST_BLAST_DATA = 'test_data/input/test_blast_PBP_alleles.out'
    TEST_BLAST_IMPERFECT_DATA = 'test_data/input/test_blast_imperfect_PBP_alleles.out'
    TEST_SEQ_DATA = 'test_data/input/test_GBS1A-1.faa'
    TEST_OUTPUT_PREFIX = 'test_data/output/GBS1A-1'
    TEST_OUTPUT_IDENTICAL_ALLELES = TEST_OUTPUT_PREFIX + '_existing_allele.txt'
    TEST_OUTPUT_NEW_ALLELES = TEST_OUTPUT_PREFIX + '_new_allele.faa'


    def test_read_blast_out(self):
        """
        Test output of blast
        """
        params_list = [
            ('.26077_6_118.11:39458-40418(+)', 0, ['1||GBS_1A', '100.000', '320', '0', '0', '1', '320', '1', '320', '0.0', '652']),
            ('.26077_6_118.11:39458-40418(+)', 1, ['138||GBS_1A', '99.688', '320', '1', '0', '1', '320', '1', '320', '0.0', '651'])
        ]
        for params in params_list:
            blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
            actual = blast_data_to_process.get_data()[params[0]][params[1]]
            self.assertEqual(actual, params[2])


    def test_get_best_blast_hit(self):
        """
        Test the best blast hit of each beta lactam
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        actual = blast_data_to_process.get_best_hit()['.26077_6_118.11:39458-40418(+)']

        self.assertEqual(actual, ['1||GBS_1A', '100.000', '320', '0', '0', '1', '320', '1', '320', '0.0', '652'])


    def test_read_seq_data(self):
        """
        Test output of sequence data
        """
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        actual = seq_data_to_analyse.get_data()

        self.assertEqual(actual['.26077_6_118.11:39458-40418(+)'],'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY')


    def test_calculate_seq_length(self):
        """
        Test the sequence lengths output
        """
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        actual = seq_data_to_analyse.calculate_seq_length()

        self.assertEqual(actual['.26077_6_118.11:39458-40418(+)'], 320)


    def test_get_identical_allele(self):
        """
        Test detection of identical alleles
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()

        self.assertEqual(get_identical_allele(best_blast_hit), ['Contig\tPBP_allele\n.26077_6_118.11:39458-40418(+)\t1||GBS_1A\n'])


    def test_get_identical_allele_with_no_identical_hits(self):
        """
        Test detection of identical alleles
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_IMPERFECT_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()

        self.assertEqual(get_identical_allele(best_blast_hit), [])


    def test_get_imperfect_allele(self):
        """
        Test detection of imperfect alleles
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_IMPERFECT_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)

        self.assertEqual(get_imperfect_allele(best_blast_hit, seq_data_to_analyse), ['>.26077_6_118.11:39458-40418(+)\nDIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n'])


    def test_get_imperfect_allele_with_identical_hits(self):
        """
        Test detection of imperfect alleles
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)

        self.assertEqual(get_imperfect_allele(best_blast_hit, seq_data_to_analyse), [])


    def test_write_content_of_identical_allele(self):
        """
        Test contents of written file
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()
        identical_allele = get_identical_allele(best_blast_hit)
        write_content(identical_allele, self.TEST_OUTPUT_IDENTICAL_ALLELES)
        fo = open(self.TEST_OUTPUT_IDENTICAL_ALLELES, 'r')

        self.assertEqual(fo.readlines(), ['Contig\tPBP_allele\n', '.26077_6_118.11:39458-40418(+)\t1||GBS_1A\n'])


    def test_write_content_of_imperfect_allele(self):
        """
        Test contents of written file
        """
        blast_data_to_process = BlastData(self.TEST_BLAST_IMPERFECT_DATA)
        best_blast_hit = blast_data_to_process.get_best_hit()
        seq_data_to_analyse = SeqData(self.TEST_SEQ_DATA)
        imperfect_allele = get_imperfect_allele(best_blast_hit, seq_data_to_analyse)
        write_content(imperfect_allele, self.TEST_OUTPUT_NEW_ALLELES)
        fo = open(self.TEST_OUTPUT_NEW_ALLELES, 'r')

        self.assertEqual(fo.readlines(), ['>.26077_6_118.11:39458-40418(+)\n', 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n'])


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--blast_out_file', 'blast_out_file', '--query_fasta', 'fasta_file',
            '--output_prefix', 'out_prefix'])
        self.assertEqual(actual, argparse.Namespace(blast_out='blast_out_file',
        fasta_qu='fasta_file', output='out_prefix'))


    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-b', 'blast_out_file', '-f', 'fasta_file', '-o', 'out_prefix'])
        self.assertEqual(actual, argparse.Namespace(blast_out='blast_out_file',
        fasta_qu='fasta_file', output='out_prefix'))


    @patch('bin.get_pbp_alleles.get_arguments')
    def test_main_with_identical_alleles(self, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.blast_out = self.TEST_BLAST_DATA
        args.fasta_qu = self.TEST_SEQ_DATA
        args.output = self.TEST_OUTPUT_PREFIX + '_MAIN'
        main()
        fo = open(self.TEST_OUTPUT_PREFIX + '_MAIN_existing_allele.txt', 'r')

        self.assertEqual(fo.readlines(), ['Contig\tPBP_allele\n', '.26077_6_118.11:39458-40418(+)\t1||GBS_1A\n'])


    @patch('bin.get_pbp_alleles.get_arguments')
    def test_main_with_imperfect_alleles(self, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        args.blast_out = self.TEST_BLAST_IMPERFECT_DATA
        args.fasta_qu = self.TEST_SEQ_DATA
        args.output = self.TEST_OUTPUT_PREFIX + '_MAIN'
        main()
        fo = open(self.TEST_OUTPUT_PREFIX + '_MAIN_new_allele.faa', 'r')

        self.assertEqual(fo.readlines(), ['>.26077_6_118.11:39458-40418(+)\n', 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n'])
