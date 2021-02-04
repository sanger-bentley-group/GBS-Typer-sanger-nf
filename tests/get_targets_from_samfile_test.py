import argparse
import unittest
from unittest.mock import patch, call

from bin.get_targets_from_samfile import get_targets, in_line, write_sam_file, write_target_sam_files, get_arguments



class TestProcessResults(unittest.TestCase):

    TEST_TARGETS = 'test_data/input/seqs_of_interest.txt'
    TEST_SAM = 'test_data/input/get_targets_test.sam'

    def test_get_targets(self):
        actual = get_targets(self.TEST_TARGETS)
        self.assertEqual(actual, ['7__PARCGBS__PARCGBS-1__7',
            '5__GYRAGBS__GYRAGBS-1__5',
            '11__23S1__23S1-1__11',
            '12__23S3__23S3-3__12',
            '16__RPOBgbs__RPOBgbs-1__16',
            '17__RPOBgbs__RPOBgbs-2__17',
            '18__RPOBgbs__RPOBgbs-3__18',
            '19__RPOBgbs__RPOBgbs-4__19'])

    def test_in_line(self):
        target='19__RPOBgbs__RPOBgbs-4__19'
        line='@HD	VN:1.0	SO:unsorted\n'
        self.assertTrue(in_line(line, target))

        line='@SQ	SN:19__RPOBgbs__RPOBgbs-4__19	LN:33\n'
        self.assertTrue(in_line(line, target))

        line='@PG	ID:bowtie2	PN:bowtie2	VN:2.2.9	CL:"/opt/bowtie2-2.2.9/bowtie2-align-s --wrapper basic-0 -q --very-sensitive-local -a -x GBS_Res_Gene-DB_Final.fasta --passthrough -1 /tmp/29.inpipe1 -2 /tmp/29.inpipe2"'
        self.assertTrue(in_line(line, target))

        line='@SQ	SN:1__CAT__CAT-1__1	LN:100\n'
        self.assertFalse(in_line(line, target))

        line='HX6_25292:2:2113:29041:9150	145	19__RPOBgbs__RPOBgbs-4__19	9	255	25M126S	18__RPOBgbs__RPOBgbs-3__18	660	\
        TCCAAAAGCACCATATGTTGGTACTGGTATGGAGTATCAAGCAGCCCACGATTCAGGTGCAGCTGTGATTGCTAAACATGACGGTCGTGTTATTTTTTCAGATGCTGAAAAAGTTGAAGTGCGTCGCGAAGATGGTTCTCTTGATGTTTAT	\
        7JJFF<FAA<<JJJJJJJJJJJFFFFJJFJJJAJFAJJJJJJJJFFFFJJJJFJJJJJJJFFJJJJJJJJJJJJJJJFFJJFFJJJJJJAJJJJJJJFFFJFJJJJJJJJJJJJJJJFJF-JJJJJJJJJJJJJJJJJJJJJJJFJFF<AA	\
        AS:i:50	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:25	YS:i:86	YT:Z:DP\n'
        self.assertTrue(in_line(line, target))

        line='HX6_25292:2:1216:15767:57389	137	2__ERMB__ERMB-1__2	1	255	80S71M	=	1	0	\
        GAAGGATTCTACAAGCGTACCTTGGATATTCACCGAACACTAGGGTTGCTCTTGCACACTCAAGTCTCGATTCAGCAATTGCTTAAGCTGCCAGCGGAATGCTTTCATCCTAAACCAAAAGTAAACAGTGTCTTAATAAAACTTACCCGCC	\
        AA<FFAJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJFJJJJJFJJFAFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJJJJJJJJJJJJJJJJFFFJJFJAFFFAJFJJAJAJ<JJ-	\
        AS:i:142	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:71	YT:Z:UP\n'
        self.assertFalse(in_line(line, target))

    def test_write_sam_file(self):
        write_sam_file(self.TEST_SAM, '12__23S3__23S3-3__12', '26189_8#5', 'test_data/output/CHECK_')
        f = open('test_data/output/CHECK_12__23S3__23S3-3__12_26189_8#5_seq.sam', "r")
        actual = "".join(f.readlines())
        self.assertEqual(actual, """@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:12__23S3__23S3-3__12\tLN:60\n@PG\tID:bowtie2\nHX4_26077:6:2110:21704:24005\t153\t12__23S3__23S3-3__12\t1\n""")


    @patch('bin.get_targets_from_samfile.write_sam_file')
    def test_write_target_sam_files(self, mock_write_sam_file):
        targets = ['7__PARCGBS__PARCGBS-1__7',
            '5__GYRAGBS__GYRAGBS-1__5',
            '11__23S1__23S1-1__11',
            '12__23S3__23S3-3__12',
            '16__RPOBgbs__RPOBgbs-1__16',
            '17__RPOBgbs__RPOBgbs-2__17',
            '18__RPOBgbs__RPOBgbs-3__18',
            '19__RPOBgbs__RPOBgbs-4__19']
        write_target_sam_files(targets, self.TEST_SAM, '26189_8#5', 'test_data/output/CHECK_')
        mock_write_sam_file.assert_has_calls([
            call(self.TEST_SAM, '7__PARCGBS__PARCGBS-1__7', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '5__GYRAGBS__GYRAGBS-1__5', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '11__23S1__23S1-1__11', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '12__23S3__23S3-3__12', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '16__RPOBgbs__RPOBgbs-1__16', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '17__RPOBgbs__RPOBgbs-2__17', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '18__RPOBgbs__RPOBgbs-3__18', '26189_8#5', 'test_data/output/CHECK_'),
            call(self.TEST_SAM, '19__RPOBgbs__RPOBgbs-4__19', '26189_8#5', 'test_data/output/CHECK_')
            ], any_order = False)

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--sam_file', 'sam_file', '--target_file', 'target_file',
            '--id', 'id', '--output_prefix', 'output'])
        self.assertEqual(actual,
                         argparse.Namespace(sam='sam_file', target='target_file', id='id', output='output'))
