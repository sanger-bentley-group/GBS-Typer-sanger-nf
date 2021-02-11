import argparse
import unittest
from unittest.mock import patch, call, ANY
from bin.get_target2MIC import MicPredictor, get_arguments, check_lists_equal, run, extract_pb2x_id


class TestGetTarget2MIC(unittest.TestCase):
    """ Test class for the get_target2MIC module """

    def setUp(self) -> None:
        self.__under_test = MicPredictor()

    def test_initialise(self) -> None:
        self.__under_test.initialise()
        self.assertEqual({}, self.__under_test.output)
        self.assertEqual(self.__under_test.PBP2X_DEFAULT, self.__under_test.pbp2x)
        self.assertEqual('neg', self.__under_test.res_dict['PBP'])
        self.assertEqual('neg', self.__under_test.res_dict['EC'])
        self.assertEqual('neg', self.__under_test.res_dict['FQ'])
        self.assertEqual('neg', self.__under_test.res_dict['OTHER'])
        self.assertEqual('neg', self.__under_test.res_dict['TET'])

    def test_check_lists_equal(self):
        self.assertTrue(check_lists_equal(['ab', 'cd', 'ef'], ['cd', 'ef', 'ab']))
        self.assertFalse(check_lists_equal(['cd', 'ef'], ['cd', 'ef', 'ab']))
        self.assertFalse(check_lists_equal(['ab', 'cd', 'ef'], ['cd', 'ab']))
        self.assertFalse(check_lists_equal(['ab', 'cd', 'ef'], ['se', 'ab', 'cd']))

    def test_get_arguments(self) -> None:
        actual = get_arguments().parse_args(
            ['--res_file', 'res_file', '--pbp_file', 'pbp_file', '--output', 'out'])
        self.assertEqual(argparse.Namespace(res_file='res_file', pbp_file='pbp_file', output_file='out'), actual)

    def test_get_arguments_short_options(self) -> None:
        actual = get_arguments().parse_args(['-r', 'res_file', '-p', 'pbp_file', '-o', 'out'])
        self.assertEqual(argparse.Namespace(
                             res_file='res_file', pbp_file='pbp_file', output_file='out'), actual)

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

    def test_extract_pb2x_id(self) -> None:
        self.assertEqual('NA', extract_pb2x_id('FOOBAR'))
        self.assertEqual('6', extract_pb2x_id('FOOBAR1:FOOBAR2:6'))
        self.assertEqual('NA', extract_pb2x_id('FOOBAR1:FOOBAR2'))
        self.assertEqual('NA', extract_pb2x_id(''))

    def test_run_prediction(self):
        res_dict = {
            'EC': '*R23S1*:*RPLD1*:foobar',
            'FQ': 'foobar:PARC-D10G',
            'TET': 'foobar1:*TET*:foobar2',
            'OTHER': 'neg',
        }
        self.__under_test.run_prediction(res_dict, self.__under_test.PBP2X_DEFAULT)

    @patch('lib.file_utils.FileUtils.read_delimited_id_file_with_hdrs')
    @patch('lib.file_utils.FileUtils.write_output')
    @patch('bin.get_target2MIC.get_arguments')
    def test_run(self, mock_get_arguments, mock_write_output, mock_read_file) -> None:
        """
        Test the main run method
        """
        mock_read_file.side_effect = [
            (
                ['EC', 'FQ', 'TET', 'OTHER'],
                {
                    '26189_8#338': {
                        'EC': '*R23S1*:*RPLD1*:foobar',
                        'FQ': 'foobar:PARC-D10G',
                        'TET': 'foobar1:*TET*:foobar2',
                        'OTHER': 'neg',
                    },
                    '26189_8#339': {
                        'EC': 'ec1',
                        'FQ': 'fq1',
                        'TET': 'foobar1:*TET*:foobar2',
                        'OTHER': 'neg',
                    }
                }
            ),
            (
                ['Contig', 'PBP_allele'],
                {
                    '26189_8#338': {
                        'Contig': '.26189_8_338.6:44447-45407(+)',
                        'PBP_allele': '1||GBS_1A',
                    },
                    '26189_8#339': {
                        'Contig': '.26189_8_338.6:44447-45407(+)',
                        'PBP_allele': '1||GBS_1A',
                    }
                }
            ),
        ]

        expected_output = \
            'ID\tPBP\tTET\tEC\tFQ\tOTHER\n' + \
            '26189_8#338\tFlag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,NA,NA,NA,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag\tfoobar1:*TET*:foobar2,>=,8,R\t*R23S1*:*RPLD1*:foobar,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag\tfoobar:PARC-D10G,NA,NA,NA,=,4,I\tneg,<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U\n' + \
            "26189_8#339\tFlag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,NA,NA,NA,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag\tfoobar1:*TET*:foobar2,>=,8,R\tec1,<=,0.25,S,<=,0.25,S,<=,2.0,S,<=,1.0,S,neg\tfq1,NA,NA,NA,Flag,Flag,Flag\tneg,<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U"

        args = mock_get_arguments.return_value.parse_args()
        args.res_file = 'res_file'
        args.pbp_file = 'pbp_file'
        args.output_file = 'out'

        run()

        mock_write_output.assert_has_calls([call(expected_output, 'out')])



