from unittest import TestCase

import os
import pandas as pd
import numpy as np

from collections import defaultdict
from bin.include_qc import *


class IncludeQCTest(TestCase):
    TEST_DATA_DIR = 'test_data'
    TEST_QC_COMPLETE = f'{TEST_DATA_DIR}/input/qc_report_complete.tab'
    TEST_QC_SUMMARY = f'{TEST_DATA_DIR}/input/qc_report_summary.tab'
    TEST_IN_SILICO = f'{TEST_DATA_DIR}/input/in_silico_data.txt'
    TEST_HEADERS = 'headers.json'
    TEST_OUTPUT = f'{TEST_DATA_DIR}/output/final_output.txt'

    header_dict = read_header_json(TEST_HEADERS)

    def test_get_content(self):
        actual = get_content(self.TEST_QC_SUMMARY, self.header_dict["id"][0], list(self.header_dict["qc_summary"].keys()), list(self.header_dict["qc_summary"].values()))

        self.assertEqual(actual.to_dict(), {
        'lane_id': {0: 'test_lane1',
                    1: 'test_lane2',
                    2: 'test_lane3',
                    3: 'test_lane4'},
        'overall_qc_status': {0: 'PASS', 1: 'FAIL', 2: 'FAIL', 3: np.nan}})

    def test_create_report(self):
        actual = create_report(self.TEST_QC_SUMMARY, self.TEST_QC_COMPLETE, self.TEST_IN_SILICO, self.header_dict)

        self.maxDiff = None
        self.assertEqual(actual.to_dict()['23S1_SNP'], {0: '*', 1: np.nan, 2: np.nan, 3: np.nan})

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--qc_summary', 'qc_summary', '--qc_complete', 'qc_complete', '--headers', 'headers', '--in_silico', 'in_silico', '--output', 'output'])
        self.assertEqual(actual, argparse.Namespace(qc_summary='qc_summary', qc_complete='qc_complete', header='headers', in_silico='in_silico', output='output'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-s', 'qc_summary', '-c', 'qc_complete', '-j', 'headers', '-i', 'in_silico', '-o', 'output'])
        self.assertEqual(actual, argparse.Namespace(qc_summary='qc_summary', qc_complete='qc_complete', header='headers', in_silico='in_silico', output='output'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--qc_summary', self.TEST_QC_SUMMARY, '--qc_complete', self.TEST_QC_COMPLETE, '--in_silico', self.TEST_IN_SILICO, '--headers', self.TEST_HEADERS, '--output', self.TEST_OUTPUT])

        main(args)

        file = open(self.TEST_OUTPUT, "r")
        actual = "".join(file.readlines())
        #os.remove(self.TEST_OUTPUT)

        self.assertEqual(actual, """lane_id	overall_qc_status	rel_abundance_status	contig_no_status	contig_gc_content_status	genome_len_status	cov_depth_status	cov_breadth_status	HET_SNPs_status	cps_type	ST	adhP	pheS	atr	glnA	sdhA	glcK	tkt	aac(6')-aph(2'')	ant(6-Ia)	aph(3'-III)	aadE	cat(pc194)	catQ	ermA	ermB	ermT	lnuB	lnuC	lsaC	lsaE	mefA	msrD	tetB	tetL	tetM	tetW	tetO	tetS	tetO32O	tetOW	tetOW32O	tetOW32OWO	tetOWO	tetSM	tetW32O	alp1	alp2/3	alpha	hvgA	PI1	PI2A1	PI2A2	PI2B	rib	srr1	srr2	23S1_SNP	23S3_SNP	gyrA_SNP	parC_SNP\ntest_lane1	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	II	1.0	1.0	1.0	2.0	1.0	1.0	2.0	2.0	neg	neg	neg	neg	neg	neg	pos	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	neg	pos	neg	neg	pos	pos	neg	neg	neg	pos	neg	*	*	*	*\ntest_lane2	FAIL	FAIL	FAIL	FAIL	FAIL	FAIL	FAIL	FAIL																																																				\ntest_lane3	FAIL	FAIL	PASS	PASS	PASS	PASS	PASS	PASS																																																				\ntest_lane4			FAIL	FAIL	FAIL	FAIL	FAIL	FAIL																																																				
""")
