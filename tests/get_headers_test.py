from unittest import TestCase

from lib.get_headers import *

class GetHeaders(TestCase):
    TEST_HEADER_FILE = 'tests/test_data/test_headers.json'

    def test_read_header_json(self):
        actual = read_header_json(self.TEST_HEADER_FILE)

        self.assertEqual(actual, {
        "contig_number": [
            "contig_no",
            "contig_no_status"
        ],
        'relative_abundance': [
            'rel_abundance',
            'rel_abundance_status'
        ],
        "contig_gc_content": [
            "gc_content",
            "gc_content_status"
        ],
        "genome_length": [
            "genome_len",
            "genome_len_status"
        ],
        "coverage_depth": [
            "cov_depth",
            "cov_depth_status"
        ],
        "coverage_breadth": [
            "cov_breadth",
            "cov_breadth_status"
        ],
        "het_snps": [
            "HET_SNPs",
            "HET_SNPs_status"
        ]})
