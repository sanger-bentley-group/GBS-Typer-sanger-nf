import argparse
import unittest
from unittest.mock import patch, call, ANY
from bin.process_surface_typer_results import get_arguments, run,  \
     derive_presence_absence, update_protein_presence_absence, featureCol, binFeatureCol, variantLookup


class TestProcessSurfaceTyperResults(unittest.TestCase):
    TEST_LANE = "26189_8#338"
    TEST_GBS_FULLGENES_RESULTS_FILE = \
        "test_data/input/" + TEST_LANE + "_SURFACE__fullgenes__GBS_Surface_Gene-DB_Final__results.txt"

    MIN_DEPTH = 30

    def setUp(self) -> None:

        self.features = {
            'ALPH': 'neg',
            'SRR': 'neg',
            'PILI': 'neg',
            'HVGA': 'neg',
        }

        self.bin_features = {
            'HVGA': '-',
            'PI1': '-',
            'PI2A1': '-',
            'PI2A2': '-',
            'PI2B': '-',
            'SRR1': '-',
            'SRR2': '-',
            'ALP1': '-',
            'ALP23': '-',
            'ALPHA': '-',
            'RIB': '-',
        }

    @patch('bin.process_surface_typer_results.update_protein_presence_absence')
    def test_derive_presence_absence(self, mock_update_protein_presence_absence):
        derive_presence_absence(
            self.TEST_GBS_FULLGENES_RESULTS_FILE, self.MIN_DEPTH, mock_update_protein_presence_absence)
        mock_update_protein_presence_absence.assert_has_calls([
            call('SRR1', 'SRR1-150', self.MIN_DEPTH, 138.384, featureCol, binFeatureCol, variantLookup),
            call('ALP23', 'ALP23-1', self.MIN_DEPTH, 166.103, featureCol, binFeatureCol, variantLookup),
            call('PI1', 'PI1-1', self.MIN_DEPTH, 137.787, featureCol, binFeatureCol, variantLookup),
            call('PI2A1', 'PI2A1-1', self.MIN_DEPTH, 126.262, featureCol, binFeatureCol, variantLookup),
            call('PI2A3', 'PI2A3-1', self.MIN_DEPTH, 100.967, featureCol, binFeatureCol, variantLookup),
            call('HVGA', 'HVGA1', self.MIN_DEPTH, 126.222, featureCol, binFeatureCol, variantLookup)])

    def test_update_protein_presence_absence_SRR(self):
        update_protein_presence_absence(
            'SRR1', 'SRR1-150', self.MIN_DEPTH, 100.0, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'neg', 'SRR': 'SRR1', 'PILI': 'neg', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '-',
                            'PI2B':  '-',
                            'SRR1':  '+',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '-',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_update_protein_presence_absence_ALPH(self):
        update_protein_presence_absence(
            'ALP23', 'ALP23-1', self.MIN_DEPTH, 100.0, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'ALP23', 'SRR': 'neg', 'PILI': 'neg', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '-',
                            'PI2B':  '-',
                            'SRR1':  '-',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '+',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_update_protein_presence_absence_PI1(self):
        update_protein_presence_absence(
            'PI1', 'PI2A2', self.MIN_DEPTH, 100.0, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'neg', 'SRR': 'neg', 'PILI': 'PI1', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '+',
                            'PI2B':  '-',
                            'SRR1':  '-',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '-',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_update_protein_presence_absence_HVGA(self):
        update_protein_presence_absence(
            'PI1', 'PI2A2', self.MIN_DEPTH, 100.0, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'neg', 'SRR': 'neg', 'PILI': 'PI1', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '+',
                            'PI2B':  '-',
                            'SRR1':  '-',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '-',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_update_protein_presence_absence_min_depth(self):
        update_protein_presence_absence(
            'PI1', 'PI2A2', self.MIN_DEPTH, 5.000, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'neg', 'SRR': 'neg', 'PILI': 'neg', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '-',
                            'PI2B':  '-',
                            'SRR1':  '-',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '-',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_update_protein_presence_absence_OTHER(self):
        update_protein_presence_absence(
            'FOOBAR', 'FOOBAR1-1', self.MIN_DEPTH, 110.000, self.features, self.bin_features, variantLookup)

        self.assertEqual(self.features, {'ALPH': 'neg', 'SRR': 'neg', 'PILI': 'neg', 'HVGA': 'neg'})
        self.assertEqual(self.bin_features,
                         {
                            'HVGA':  '-',
                            'PI1':   '-',
                            'PI2A1': '-',
                            'PI2A2': '-',
                            'PI2B':  '-',
                            'SRR1':  '-',
                            'SRR2':  '-',
                            'ALP1':  '-',
                            'ALP23': '-',
                            'ALPHA': '-',
                            'RIB':   '-'})

    def test_get_arguments(self):
        actual = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_output_name', '--surface_db', 'surface_db', '--output', 'outfile',
             '--min_read_depth', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(
                            fullgenes_file_id='srst2_output_name', db='surface_db', output='outfile', min_depth=30.0))

    def test_get_arguments_short_options(self):
        actual = get_arguments().parse_args(['-s', 'srst2_output_name', '-b', 'sero_db', '-o', 'outfile', '-d', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(
                             fullgenes_file_id='srst2_output_name', db='sero_db', output='outfile', min_depth=30.0))

    @patch('bin.process_surface_typer_results.derive_presence_absence')
    @patch('lib.file_utils.FileUtils.create_output_contents')
    @patch('lib.file_utils.FileUtils.write_output')
    def test_run(self, mock_write_output, mock_create_output_contents, mock_derive_presence_absence):

        # Given

        args = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_gbs_out',
             '--surface_db', 'gbs_surface_db.fasta',
             '--min_read_depth', '40.0', '--output_prefix', 'output'])
        mock_create_output_contents.return_value = 'foobar'

        # When

        run(args)

        # Then

        self.assertEqual(
            mock_derive_presence_absence.call_args_list,
            [call('srst2_gbs_out__fullgenes__gbs_surface_db__results.txt', 40.0, ANY)])
        mock_create_output_contents.assert_has_calls([call(ANY), call(ANY)])
        mock_write_output.assert_has_calls([
            call(ANY, args.output + '_surface_protein_variants_sample.txt'),
            call(ANY, args.output + '_surface_protein_incidence_sample.txt')
        ], any_order=False)
