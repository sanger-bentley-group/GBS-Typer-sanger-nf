import argparse
import unittest
from unittest.mock import patch, call, ANY
from bin.get_target2MIC import MicPredictor


class TestGetTarget2MIC(unittest.TestCase):
    """ Test class for the get_target2MIC module """

    TEST_LANE = "26189_8#338"
    TEST_GBS_FULLGENES_RESULTS_FILE = \
        "test_data/" + TEST_LANE + "_SURFACE__fullgenes__GBS_Surface_Gene-DB_Final__results.txt"

    def setUp(self) -> None:
        self.__under_test = MicPredictor()

    def test_initialise(self) -> None:
        self.__under_test.initialise()
        self.assertEqual({}, self.__under_test.output)
        self.assertEqual(0.0, self.__under_test.pbp2x)
        self.assertEqual('neg', self.__under_test.res_dict['EC'])
        self.assertEqual('neg', self.__under_test.res_dict['FQ'])
        self.assertEqual('neg', self.__under_test.res_dict['OTHER'])
        self.assertEqual('neg', self.__under_test.res_dict['TET'])

    def test_check_lists_equal(self):
        self.assertTrue(self.__under_test.check_lists_equal(['ab', 'cd', 'ef'], ['cd', 'ef', 'ab']))
        self.assertFalse(self.__under_test.check_lists_equal(['cd', 'ef'], ['cd', 'ef', 'ab']))
        self.assertFalse(self.__under_test.check_lists_equal(['ab', 'cd', 'ef'], ['cd', 'ab']))
        self.assertFalse(self.__under_test.check_lists_equal(['ab', 'cd', 'ef'], ['se', 'ab', 'cd']))

    def test_get_arguments(self) -> None:
        actual = self.__under_test.get_arguments().parse_args(
            ['--res_file', 'res_file', '--output_prefix', 'out'])
        self.assertEqual(argparse.Namespace(res_file='res_file', output='out'), actual)

    def test_get_arguments_short_options(self) -> None:
        actual = self.__under_test.get_arguments().parse_args(['-r', 'res_file', '-o', 'out'])
        self.assertEqual(argparse.Namespace(
                             res_file='res_file', output='out'), actual)

    def test_update_pbp_category_small_pbp2x(self) -> None:
        self.__under_test.pbp2x = 1
        self.__under_test.update_pbp_category()

        self.assertEqual('<=,0.5,U,<=,8.0,U,<=,0.12,S,<=,0.12,S,<=,0.12,S,NA,NA,NA,<=,0.25,S,<=,0.12,S,<=,0.12,S',
                         self.__under_test.output['PBP'])

    def test_update_pbp_category_large_pbp2x(self) -> None:
        self.__under_test.pbp2x = 10
        self.__under_test.update_pbp_category()

        self.assertEqual(
            'Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,NA,NA,NA,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['PBP']
        )

    def test_update_er_cl_category_neg(self) -> None:
        self.__under_test.res_dict['EC'] = 'neg'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            'neg,<=,0.25,S,<=,0.25,S,<=,2.0,S,<=,1.0,S,neg',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg(self) -> None:
        self.__under_test.res_dict['EC'] = '**R23S1**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**R23S1**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_r23sd1(self) -> None:
        self.__under_test.res_dict['EC'] = '**R23S1**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**R23S1**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rpld1(self) -> None:
        self.__under_test.res_dict['EC'] = '**RPLD1**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**RPLD1**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rplv(self) -> None:
        self.__under_test.res_dict['EC'] = '**RPLV**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**RPLV**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rpld1(self) -> None:
        self.__under_test.res_dict['EC'] = '**RPLD1**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**RPLD1**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rpld1(self) -> None:
        self.__under_test.res_dict['EC'] = '**rplv**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**rplv**,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_mef(self) -> None:
        self.__under_test.res_dict['EC'] = '**MEF**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**MEF**,>=,1,R,<=,0.25,S,<=,2.0,S,<=,1.0,S,neg',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_lsa(self) -> None:
        self.__under_test.res_dict['EC'] = '**LSA**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**LSA**,<=,0.25,S,>=,1,R,<=,2.0,S,<=,1.0,S,neg',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_lnu(self) -> None:
        self.__under_test.res_dict['EC'] = '**LNU**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**LNU**,<=,0.25,S,Flag,Flag,Flag,<=,2.0,S,<=,1.0,S,neg',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_erm(self) -> None:
        self.__under_test.res_dict['EC'] = '**erm**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**erm**,>=,1,R,>=,1,R,<=,2.0,S,<=,1.0,S,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_mef_and_lsa(self) -> None:
        self.__under_test.res_dict['EC'] = '**mef**:**lsa**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**mef**:**lsa**,>=,1,R,>=,1,R,<=,2.0,S,<=,1.0,S,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_mef_and_lnu(self) -> None:
        self.__under_test.res_dict['EC'] = '**mef**:**LNU**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '**mef**:**LNU**,>=,1,R,Flag,Flag,Flag,<=,2.0,S,<=,1.0,S,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_erm_and_lsa(self) -> None:
        self.__under_test.res_dict['EC'] = '$$BAR$$:**ERM**:**LSA**:**FOO%%'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '$$BAR$$:**ERM**:**LSA**:**FOO%%,>=,1,R,>=,1,R,<=,2.0,S,Flag,Flag,Flag,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_erm_and_lnu(self) -> None:
        self.__under_test.res_dict['EC'] = '$$BAR$$:**ERM**:**LNU**:**FOO%%'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            '$$BAR$$:**ERM**:**LNU**:**FOO%%,>=,1,R,>=,1,R,<=,2.0,S,<=,1.0,S,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rpld1_erm_and_lnu(self) -> None:
        self.__under_test.res_dict['EC'] = 'RPLD1:ERM:LNU**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            'RPLD1:ERM:LNU**,>=,1,R,>=,1,R,Flag,Flag,Flag,<=,1.0,S,pos',
            self.__under_test.output['EC']
        )

    def test_update_er_cl_category_non_neg_rpld1_mef(self) -> None:
        self.__under_test.res_dict['EC'] = 'RPLD1:FOO:MEF1**'
        self.__under_test.update_er_cl_category()

        self.assertEqual(
            'RPLD1:FOO:MEF1**,>=,1,R,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['EC']
        )

    def test_update_tet_category_neg(self) -> None:
        self.__under_test.res_dict['TET'] = 'neg'
        self.__under_test.update_tet_category()
        self.assertEqual('neg,<=,2.0,S', self.__under_test.output['TET'])

    def test_update_tet_category_non_neg1(self) -> None:
        self.__under_test.res_dict['TET'] = 'foobar0:**TET**:foobar1:foobar2'
        self.__under_test.update_tet_category()
        self.assertEqual('foobar0:**TET**:foobar1:foobar2,>=,8,R', self.__under_test.output['TET'])

    def test_update_tet_category_non_neg2(self) -> None:
        self.__under_test.res_dict['TET'] = 'foobar0:foobar1:foobar2'
        self.__under_test.update_tet_category()
        self.assertEqual('foobar0:foobar1:foobar2,neg', self.__under_test.output['TET'])

    def test_update_gyra_parc_category_neg(self) -> None:
        self.__under_test.res_dict['FQ'] = 'neg'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('neg,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

    def test_update_gyra_parc_category_non_neg_gyra1(self) -> None:
        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-S6F**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-S6F**,NA,NA,NA,>=,8,R', self.__under_test.output['FQ'])

    def test_update_gyra_parc_category_non_neg_gyra2(self) -> None:
        self.__under_test.res_dict['FQ'] = 'FOO:**FOOBAR**:**PARC-S6Y**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**FOOBAR**:**PARC-S6Y**,NA,NA,NA,=,4,I', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-D10G**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-D10G**,NA,NA,NA,=,4,I', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-D10Y**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-D10Y**,NA,NA,NA,=,4,I', self.__under_test.output['FQ'])

    def test_update_gyra_parc_category_non_neg_gyra3(self) -> None:
        self.__under_test.res_dict['FQ'] = 'FOO:**FOOBAR**:**PARC-D10A**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**FOOBAR**:**PARC-D10A**,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-D10N**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-D10N**,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-D5N**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-D5N**,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**PARC-S6F**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**PARC-S6F**,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-S7P**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-S7P**,NA,NA,NA,<=,2,S', self.__under_test.output['FQ'])

    def test_update_gyra_parc_category_non_neg_default(self) -> None:
        self.__under_test.res_dict['FQ'] = 'FOO:**GYRA-S11L**:**PARC-S6ZX**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**GYRA-S11L**:**PARC-S6ZX**,NA,NA,NA,Flag,Flag,Flag', self.__under_test.output['FQ'])

        self.__under_test.res_dict['FQ'] = 'FOO:**PARC-D10**:**PARC-S6ZX**'
        self.__under_test.update_gyra_parc_category()
        self.assertEqual('FOO:**PARC-D10**:**PARC-S6ZX**,NA,NA,NA,Flag,Flag,Flag', self.__under_test.output['FQ'])

    def test_update_other_category_neg(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'neg'
        self.__under_test.update_other_category()
        self.assertEqual('neg,<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U', self.__under_test.output['OTHER'])

    def test_update_other_category_non_neg(self) -> None:
        self.__under_test.res_dict['OTHER'] = "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly"
        self.__under_test.update_other_category()
        self.assertEqual(
            "ant(6)-Ia:Ant6-Ia_AGly:aph(3')-III:Aph3-III_AGly:Sat4A_Agly," +
            '<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

        self.__under_test.res_dict['OTHER'] = "aph(3')-III:Aph3-III_AGly:Sat4A_AGly"
        self.__under_test.update_other_category()
        self.assertEqual(
            "aph(3')-III:Aph3-III_AGly:Sat4A_AGly," +
            '<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

        self.__under_test.res_dict['OTHER'] = "msr(D):MsrD_MLS"
        self.__under_test.update_other_category()
        self.assertEqual(
            "msr(D):MsrD_MLS," +
            '<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_cat(self) -> None:
        self.__under_test.res_dict['OTHER'] = '**CAT**:FOOBAR'
        self.__under_test.update_other_category()
        self.assertEqual(
            '**CAT**:FOOBAR,<=,1,S,<=,1,S,<=,1,U,>=,16,R,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_fola(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'FOOBAR:**FOLA**'
        self.__under_test.update_other_category()
        self.assertEqual(
            'FOOBAR:**FOLA**,<=,1,S,<=,1,S,<=,1,U,<=,4,S,Flag,Flag,Flag',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_folp(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'FOOBAR:**FOLP**'
        self.__under_test.update_other_category()
        self.assertEqual(
            'FOOBAR:**FOLP**,<=,1,S,<=,1,S,<=,1,U,<=,4,S,Flag,Flag,Flag',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_rpob(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'FOOBAR:**rpob**'
        self.__under_test.update_other_category()
        self.assertEqual(
            'FOOBAR:**rpob**,<=,1,S,<=,1,S,Flag,Flag,Flag,<=,4,S,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_van(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'FOOBAR:**VAN**'
        self.__under_test.update_other_category()
        self.assertEqual(
            'FOOBAR:**VAN**,<=,1,S,>=,2,U,<=,1,U,<=,4,S,<=,0.5,U',
            self.__under_test.output['OTHER']
        )

    def test_update_other_category_non_neg_flags(self) -> None:
        self.__under_test.res_dict['OTHER'] = 'FOOBAR1:FOOBAR2:FOOBAR3'
        self.__under_test.update_other_category()
        self.assertEqual(
            'FOOBAR1:FOOBAR2:FOOBAR3,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag',
            self.__under_test.output['OTHER']
        )

    def test_read_res_alleles_file(self) -> None:
        self.__under_test.read_res_alleles_file('\t', 'test_data/input/25292_2#85_res_alleles.txt')
        print(self.__under_test.res_dict)
        self.assertEqual(4, len(self.__under_test.res_dict))
        self.assertEqual("ERMB(ERMB-1):MEF(MEF-1):23S1:23S3", self.__under_test.res_dict['EC'])
        self.assertEqual("PARC:GYRA-V1A,M2Q,G3W,K4W", self.__under_test.res_dict['FQ'])
        self.assertEqual(
            "MsrD_MLS(MsrD_296):Ant6-Ia_AGly(Ant6-Ia_1633):Sat4A_AGly(Sat4A_586):Aph3-III_AGly(Aph3-III_268):msr(D)(msr(D)_2):ant(6)-Ia(ant(6)-Ia):aph(3')-III(aph(3')-III)",
            self.__under_test.res_dict['OTHER']
        )
        self.assertEqual("TETO(TETO-1)", self.__under_test.res_dict['TET'])

    def test_read_res_alleles_file_unexpected_format(self) -> None:
        self.assertRaises(
            ValueError,
            self.__under_test.read_res_alleles_file,
            '\t',
            'test_data/input/25292_2#85_res_alleles-badformat.txt')

    def test_read_res_alleles_file_bad_file(self) -> None:
        self.assertRaises(
            IOError,
            self.__under_test.read_res_alleles_file,
            '\t',
            'test_data/input/non-existent.txt')

    def test_search_list(self) -> None:
        test_arr = ['fred1', 'bob3', 'bob4', 'jane345', 'mary-1234', 'david_22334']

        self.assertFalse(self.__under_test.search_list(test_arr, r"foobar"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"bob"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"jim|bob|susan"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"jim|foobar|bob[3|7]"))
        self.assertFalse(self.__under_test.search_list(test_arr, r"jim|foobar|bob[1|2]"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"david_22334"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"fred1"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"fred[1|3|6]"))
        self.assertTrue(self.__under_test.search_list(test_arr, r"FRED[1|3|6]"))
        self.assertFalse(self.__under_test.search_list([], r"jane345"))

    @patch('bin.get_target2MIC.MicPredictor.read_res_alleles_file')
    @patch('lib.file_utils.FileUtils.create_output_contents')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.get_target2MIC.MicPredictor.get_arguments')
    def test_run(self, mock_get_arguments, mock_write_output, mock_create_output_contents, mock_res_file) -> None:
        res_dict = {'EC': '*R23S1*:*RPLD1*:foobar'}
        self.__under_test.pbp2x = 2

        args = self.__under_test.get_arguments().parse_args(['--res_file', 'res_file', '--output_prefix', 'out'])

        mock_res_file.return_value = res_dict
        self.__under_test.run()

        mock_create_output_contents.assert_has_calls([call(ANY)])
        mock_write_output.assert_has_calls([call(ANY, args.output + self.__under_test.MIC_PREDICTIONS_OUTPUT_FILE)])



