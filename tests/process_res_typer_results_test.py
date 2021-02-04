import argparse
import unittest
from unittest.mock import patch, call, ANY

from bin.process_res_typer_results import get_arguments, codon2aa, derive_presence_absence_targets, \
    derive_presence_absence_targets_for_arg_res, six_frame_translate, find_mismatches, update_presence_absence_target, \
    update_presence_absence_target_for_arg_res, drugRes_Col, get_seq_diffs, update_GBS_Res_var, update_drug_res_col_dict, \
    get_gene_names_from_consensus, get_variants, run, main, get_seq_content, \
    geneToRef, GBS_Res_var, Res_Targets, geneToClass, extract_frame_aa, EOL_SEP


class TestProcessResTyperResults(unittest.TestCase):

    TEST_LANE = "26189_8#5"
    TEST_GBS_FULLGENES_RESULTS_FILE = "test_data/input/RES_" + TEST_LANE + "__fullgenes__GBS_Res_Gene-DB_Final__results.txt"
    TEST_ARGANNOT_FULLGENES_RESULTS_FILE = "test_data/input/ARG_" + TEST_LANE + "__fullgenes__ARG-ANNOT__results.txt"
    TEST_RESFINDER_FULLGENES_RESULTS_FILE = "test_data/input/RESFI_" + TEST_LANE + "__fullgenes__ResFinder__results.txt"
    TEST_FASTA_FILE = "test_data/input/test-db.fasta"
    TEST_CONSENSUS_SEQ_FILE = "test_data/input/" + TEST_LANE + "_consensus_seq.fna"
    TEST_OUTPUT = "test_data/output/" + TEST_LANE + "_output.txt"

    MIN_DEPTH = 30

    def test_codon2aa(self):
        self.assertEqual('S', codon2aa('tca'))
        self.assertEqual('S', codon2aa('TCA'))
        self.assertEqual('S', codon2aa('TCC'))
        self.assertEqual('S', codon2aa('TCG'))
        self.assertEqual('S', codon2aa('TCT'))
        self.assertEqual('F', codon2aa('TTC'))
        self.assertEqual('F', codon2aa('TTT'))
        self.assertEqual('L', codon2aa('TTA'))
        self.assertEqual('L', codon2aa('TTG'))
        self.assertEqual('Y', codon2aa('TAC'))
        self.assertEqual('Y', codon2aa('TAT'))
        self.assertEqual('*', codon2aa('TAA'))
        self.assertEqual('*', codon2aa('TAG'))
        self.assertEqual('C', codon2aa('TGC'))
        self.assertEqual('C', codon2aa('TGT'))
        self.assertEqual('*', codon2aa('TGA'))
        self.assertEqual('W', codon2aa('TGG'))
        self.assertEqual('L', codon2aa('CTA'))
        self.assertEqual('L', codon2aa('CTC'))
        self.assertEqual('L', codon2aa('CTG'))
        self.assertEqual('L', codon2aa('CTT'))
        self.assertEqual('P', codon2aa('CCA'))
        self.assertEqual('P', codon2aa('CCC'))
        self.assertEqual('P', codon2aa('CCG'))
        self.assertEqual('P', codon2aa('CCT'))
        self.assertEqual('H', codon2aa('CAC'))
        self.assertEqual('H', codon2aa('CAT'))
        self.assertEqual('Q', codon2aa('CAA'))
        self.assertEqual('Q', codon2aa('CAG'))
        self.assertEqual('R', codon2aa('CGA'))
        self.assertEqual('R', codon2aa('CGC'))
        self.assertEqual('R', codon2aa('CGG'))
        self.assertEqual('R', codon2aa('CGT'))
        self.assertEqual('I', codon2aa('ATA'))
        self.assertEqual('I', codon2aa('ATC'))
        self.assertEqual('I', codon2aa('ATT'))
        self.assertEqual('M', codon2aa('ATG'))
        self.assertEqual('T', codon2aa('ACA'))
        self.assertEqual('T', codon2aa('ACC'))
        self.assertEqual('T', codon2aa('ACG'))
        self.assertEqual('T', codon2aa('ACT'))
        self.assertEqual('N', codon2aa('AAC'))
        self.assertEqual('N', codon2aa('AAT'))
        self.assertEqual('K', codon2aa('AAA'))
        self.assertEqual('K', codon2aa('AAG'))
        self.assertEqual('S', codon2aa('AGC'))
        self.assertEqual('S', codon2aa('AGT'))
        self.assertEqual('R', codon2aa('AGA'))
        self.assertEqual('R', codon2aa('AGG'))
        self.assertEqual('V', codon2aa('GTA'))
        self.assertEqual('V', codon2aa('GTC'))
        self.assertEqual('V', codon2aa('GTG'))
        self.assertEqual('V', codon2aa('GTT'))
        self.assertEqual('A', codon2aa('GCA'))
        self.assertEqual('A', codon2aa('GCC'))
        self.assertEqual('A', codon2aa('GCG'))
        self.assertEqual('A', codon2aa('GCT'))
        self.assertEqual('D', codon2aa('GAC'))
        self.assertEqual('D', codon2aa('GAT'))
        self.assertEqual('E', codon2aa('GAA'))
        self.assertEqual('E', codon2aa('GAG'))
        self.assertEqual('G', codon2aa('GGA'))
        self.assertEqual('G', codon2aa('GGC'))
        self.assertEqual('G', codon2aa('GGG'))
        self.assertEqual('G', codon2aa('GGT'))

        self.assertEqual('A', codon2aa('GC?'))
        self.assertEqual('G', codon2aa('GG?'))
        self.assertEqual('P', codon2aa('CC?'))
        self.assertEqual('T', codon2aa('AC?'))
        self.assertEqual('V', codon2aa('GT?'))
        self.assertEqual('R', codon2aa('CG?'))
        self.assertEqual('S', codon2aa('TC?'))

    def test_codon2aa_unknown_codon(self):

        self.assertEqual('x', codon2aa('+++'))

    def test_extract_frame_aa(self):
        self.assertEqual("IHMVIHL", extract_frame_aa("ATCCACATGGTGATTCATCTA", 1))
        self.assertEqual("IHMVIHL", extract_frame_aa("atccacatggtgattcatcta", 1))
        self.assertEqual("STW*FI", extract_frame_aa("ATCCACATGGTGATTCATCTA", 2))
        self.assertEqual("PHGDSS", extract_frame_aa("ATCCACATGGTGATTCATCTA", 3))
        self.assertEqual("*MNHHVD", extract_frame_aa("ATCCACATGGTGATTCATCTA", 4))
        self.assertEqual("*MNHHVD", extract_frame_aa("atccacatggtgattcatcta", 4))
        self.assertEqual("R*ITMW", extract_frame_aa("ATCCACATGGTGATTCATCTA", 5))
        self.assertEqual("DESPCG", extract_frame_aa("ATCCACATGGTGATTCATCTA", 6))

    def test_six_frame_translate__frame1(self):
        self.assertEqual(
            "VMGKYHPHGDSSIYEAMVRMAQWW",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGG",
                1
            )
        )

        self.assertEqual(
            "HPHGDSSIYDAMVRMSQ",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                1
            )
        )

    def test_six_frame_translate__frame2(self):
        self.assertEqual(
            "LWVNTIHMVIHLFTKQWCVWHNG",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGG",
                2
            )
        )

        self.assertEqual(
            "ILMGIPLSMTRWFVCL",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                2
            )
        )

    def test_six_frame_translate__frame3(self):
        self.assertEqual(
            "YG*IPSTW*FIYLRSNGAYGTMVEKK",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                3
            )
        )

        self.assertEqual(
            "SSWGFLYL*RDGSYVS",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                3
            )
        )

    def test_six_frame_translate__frame4(self):
        self.assertEqual(
            "FFFSTIVPYAPLLRK*MNHHVDGIYP*",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                4
            )
        )

        self.assertEqual(
            "LRHTNHRVIDRGIPMRI",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                4
            )
        )

    def test_six_frame_translate__frame5(self):
        self.assertEqual(
            "FFFPPLCHTHHCFVNR*ITMWMVFTHK",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                5
            )
        )

        self.assertEqual(
            "*DIRTIAS*IEESP*G",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                5
            )
        )

    def test_six_frame_translate__frame6(self):
        self.assertEqual(
            "FFFHHCAIRTIAS*IDESPCGWYLPI",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                6
            )
        )

        self.assertEqual(
            "ETYEPSRHR*RNPHED",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                6
            )
        )

    def test_six_frame_translate__frame_range_error(self):
        with self.assertRaises(IndexError):
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                7
            )

        with self.assertRaises(IndexError):
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                0
            )

    def test_update_presence_absence_target(self):
        depth = self.MIN_DEPTH+1


        # ============== Test misc ==================
        misc_list = ["PARC", "GYRA", "23S1", "23S3", "RPOBGBS-1"]
        for allele in misc_list:
            gbs_res_target_dict = {allele: "neg"}
            update_presence_absence_target("GENE1", "***"+allele+"***", depth, gbs_res_target_dict)
            self.assertEqual({allele: "pos"}, gbs_res_target_dict)

        # ============== Test RPOBgbs-N ==================
        gbs_res_target_dict = {"RPOBGBS-2": "neg"}
        update_presence_absence_target("GENE1", "***RPOBGBS-2***", depth, gbs_res_target_dict)
        self.assertEqual({"RPOBGBS-2": "pos"}, gbs_res_target_dict)

        # ============== Test depth ==================
        gbs_res_target_dict = {}
        update_presence_absence_target("GENE1", "***RPOBGBS-1***", self.MIN_DEPTH-1, gbs_res_target_dict)
        self.assertEqual({}, gbs_res_target_dict)

        # TODO there is a suspected bug in this perl code - see Python module

    def test_update_presence_absence_target_for_arg_res(self):
        depth = self.MIN_DEPTH+1

        # ============== Test ERMB ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERMB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***ERMB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***ERMB***]"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***ERMB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***ERMB***]:GENE2[***ERMB***]"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "pos"}, res_target_dict)

        # Check low depth
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERMB": "neg"}
        update_presence_absence_target_for_arg_res("GENE2", "***ERMB***", self.MIN_DEPTH-1, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "neg"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "neg"}, res_target_dict)

        # ============== Test TETM ==================
        drug_res_col_dict = {"TET": "neg"}
        res_target_dict = {"TETM": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***TETM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "GENE1[***TETM***]"}, drug_res_col_dict)
        self.assertEqual({"TETM": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***TETM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "GENE1[***TETM***]:GENE2[***TETM***]"}, drug_res_col_dict)
        self.assertEqual({"TETM": "pos"}, res_target_dict)

        # ============== Test CAT ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***CAT***]"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***CAT***]:GENE2[***CAT***]"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        # ============== Test LNUB ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LNUB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***LNUB***]"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***LNUB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***LNUB***]:GENE2[***LNUB***]"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        # ============== Test LSAC ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LSAC": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSAC***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***LSAC***]"}, drug_res_col_dict)
        self.assertEqual({"LSAC": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***LSAC***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***LSAC***]:GENE2[***LSAC***]"}, drug_res_col_dict)
        self.assertEqual({"LSAC": "pos"}, res_target_dict)

        # ============== Test MEFA ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"MEFA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***MEFA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***MEFA***]"}, drug_res_col_dict)
        self.assertEqual({"MEFA": "pos"}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***MEFA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***MEFA***]:GENE2[***MEFA***]"}, drug_res_col_dict)
        self.assertEqual({"MEFA": "pos"}, res_target_dict)

        # ============== Test FOSA ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"FOSA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***fosA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***fosA***]"}, drug_res_col_dict)
        self.assertEqual({"FOSA": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "CAT[***CAT***]"}
        res_target_dict = {"FOSA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***fosA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "CAT[***CAT***]:GENE1[***fosA***]"}, drug_res_col_dict)
        self.assertEqual({"FOSA": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "FOSA[***allele***]"}
        res_target_dict = {"FOSA": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***fosA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "FOSA[***allele***]:GENE1[***fosA***]"}, drug_res_col_dict)
        self.assertEqual({"FOSA": "pos"}, res_target_dict)

        # ============== Test ERMB ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERMB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***erm(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***erm(B)***]"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "LNUB[***allele***]"}
        res_target_dict = {"ERMB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***erm(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "LNUB[***allele***]:GENE1[***erm(B)***]"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERMB[***allele***]"}
        res_target_dict = {"ERMB": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***erm(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERMB[***allele***]:GENE1[***erm(B)***]"}, drug_res_col_dict)
        self.assertEqual({"ERMB": "pos"}, res_target_dict)

        # ============== Test LNUB ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***lnu(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***lnu(B)***]"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM[***allele***]"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***lnu(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM[***allele***]:GENE1[***lnu(B)***]"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "LNUB[***allele***]"}
        res_target_dict = {"LNUB": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***lnu(B)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "LNUB[***allele***]:GENE1[***lnu(B)***]"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        # ============== Test LSAC ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LSAC": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***lsa(C)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***lsa(C)***]"}, drug_res_col_dict)
        self.assertEqual({"LSAC": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM[***allele***]"}
        res_target_dict = {"LSAC": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSAC***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM[***allele***]:GENE1[***LSAC***]"}, drug_res_col_dict)
        self.assertEqual({"LSAC": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM[***allele***]"}
        res_target_dict = {"LSAC": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSAC***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM[***allele***]:GENE1[***LSAC***]"}, drug_res_col_dict)
        self.assertEqual({"LSAC": "pos"}, res_target_dict)

        # ============== Test MEFA ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"MEFA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***mef(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***mef(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MEFA": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM[***allele***]"}
        res_target_dict = {"MEFA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***mef(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM[***allele***]:GENE1[***mef(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MEFA": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM[***allele***]"}
        res_target_dict = {"MEFA": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***mef(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM[***allele***]:GENE1[***mef(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MEFA": "pos"}, res_target_dict)

        # ============== Test MPHC ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"MPHC": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***mph(C)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1[***mph(C)***]"}, drug_res_col_dict)
        self.assertEqual({"MPHC": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERMB[***allele***]"}
        res_target_dict = {"MPHC": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***mph(C)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERMB[***allele***]:GENE1[***mph(C)***]"}, drug_res_col_dict)
        self.assertEqual({"MPHC": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "MPHC[***allele***]"}
        res_target_dict = {"MPHC": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***mph(C)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "MPHC[***allele***]:GENE1[***mph(C)***]"}, drug_res_col_dict)
        self.assertEqual({"MPHC": "pos"}, res_target_dict)

        # ============== Test MSRA ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"MSRA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***msr(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRA": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "FOSA[***allele***]"}
        res_target_dict = {"MSRA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "FOSA[***allele***]:GENE1[***msr(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRA": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "MSRA[***allele***]"}
        res_target_dict = {"MSRA": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(A)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "MSRA[***allele***]:GENE1[***msr(A)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRA": "pos"}, res_target_dict)

        # ============== Test MSRD ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"MSRD": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(D)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***msr(D)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRD": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "FOSA[***allele***]"}
        res_target_dict = {"MSRD": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(D)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "FOSA[***allele***]:GENE1[***msr(D)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRD": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "MSRD[***allele***]"}
        res_target_dict = {"MSRD": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***msr(D)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "MSRD[***allele***]:GENE1[***msr(D)***]"}, drug_res_col_dict)
        self.assertEqual({"MSRD": "pos"}, res_target_dict)

        # ============== Test SUL2 ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"SUL2": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***sul2***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***sul2***]"}, drug_res_col_dict)
        self.assertEqual({"SUL2": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "FOSA[***allele***]"}
        res_target_dict = {"SUL2": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***sul2***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "FOSA[***allele***]:GENE1[***sul2***]"}, drug_res_col_dict)
        self.assertEqual({"SUL2": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "SUL2[***allele***]"}
        res_target_dict = {"SUL2": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***sul2***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "SUL2[***allele***]:GENE1[***sul2***]"}, drug_res_col_dict)
        self.assertEqual({"SUL2": "pos"}, res_target_dict)

        # ============== Test TETM ==================
        drug_res_col_dict = {"TET": "neg"}
        res_target_dict = {"TETM": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***tet(M)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "GENE1[***tet(M)***]"}, drug_res_col_dict)
        self.assertEqual({"TETM": "pos"}, res_target_dict)

        drug_res_col_dict = {"TET": "ERM[***allele***]"}
        res_target_dict = {"TETM": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***tet(M)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "ERM[***allele***]:GENE1[***tet(M)***]"}, drug_res_col_dict)
        self.assertEqual({"TETM": "pos"}, res_target_dict)

        drug_res_col_dict = {"TET": "ERM[***allele***]"}
        res_target_dict = {"TETM": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***tet(M)***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "ERM[***allele***]:GENE1[***tet(M)***]"}, drug_res_col_dict)
        self.assertEqual({"TETM": "pos"}, res_target_dict)

        # ============== Test CAT ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***CAT***]"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "ERM[***allele***]"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "ERM[***allele***]:GENE1[***CAT***]"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "CAT[***allele***]"}
        res_target_dict = {"CAT": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "CAT[***allele***]:GENE1[***CAT***]"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        # ============== Test OTHER ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {}
        update_presence_absence_target_for_arg_res("GENE1", "***FOO***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***FOO***]"}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***FOO***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1[***FOO***]:GENE2[***FOO***]"}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)

        # ============== Test depth ==================
        drug_res_col_dict = {}
        res_target_dict = {}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", self.MIN_DEPTH - 1, drug_res_col_dict, res_target_dict)
        self.assertEqual({}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)

    @patch('bin.process_res_typer_results.update_presence_absence_target')
    def test_derive_presence_absence_targets(self, mock):

        calls = [call("23S1", "23S1-1", 1135.571, ANY), call("23S3", "23S3-3", 1265.721, ANY)]
        derive_presence_absence_targets(self.TEST_GBS_FULLGENES_RESULTS_FILE)
        mock.assert_has_calls(calls, any_order=False)

    @patch('bin.process_res_typer_results.update_presence_absence_target_for_arg_res')
    def derive_presence_absence_targets_for_arg_res(self, mock):
        calls = [
            call("tet(M)", "tet(M)_12", 132.04, ANY, ANY),
            call("tet(M)", "tet(M)_4", 185.331, ANY, ANY),
            call("tet(M)", "tet(M)_10", 120.412, ANY, ANY),
        ]
        derive_presence_absence_targets_for_arg_res(self.TEST_RESFINDER_FULLGENES_RESULTS_FILE)
        mock.assert_has_calls(calls, any_order=False)

    def test_find_amino_acid_mismatches(self):
        actual = find_mismatches([], 'HPHGDSSIYDAMVRMSS', geneToRef['PARC'])
        self.assertEqual(actual, ['Q17S'])

        actual = find_mismatches([], 'HHHGDSSIYDAMVRMSS', geneToRef['PARC'])
        self.assertEqual(actual, ['P2H', 'Q17S'])

        actual = find_mismatches([], 'MMGKYHPHGDSSIYEAMVRMAQWW', geneToRef['GYRA'])
        self.assertEqual(actual, ['V1M'])

        actual = find_mismatches([], 'GGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL', geneToRef['RPOBGBS-1'])
        self.assertEqual(actual, ['F1G'])

        actual = find_mismatches([], 'SSQLVRSPGV', geneToRef['RPOBGBS-2'])
        self.assertEqual(actual, ['V1S'])

        actual = find_mismatches([], 'TTVAQANSKLNEDGTFAEEIVMGRHQGNNQEFPSSI', geneToRef['RPOBGBS-3'])
        self.assertEqual(actual, ['F1T'])

        actual = find_mismatches([], 'IIDPKAPYVGT', geneToRef['RPOBGBS-4'])
        self.assertEqual(actual, ['L1I'])

    def test_find_nucleotide_mismatches(self):
        actual = find_mismatches([], 'ATTACCCGCGACAGGACGGAAAGACCCCATGGAG', geneToRef['23S1'])
        self.assertEqual(actual, ['G1A'])

        actual = find_mismatches([], 'ATTACCCGCGACAGGACGGAAAGACCCCATGGAT', geneToRef['23S1'])
        self.assertEqual(actual, ['G1A', 'G34T'])

        actual = find_mismatches([], 'GGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG', geneToRef['23S3'])
        self.assertEqual(actual, ['C1G'])

    @patch('bin.process_res_typer_results.six_frame_translate')
    def test_get_seq_diffs(self, mock_six_frame_translate):
        mock_six_frame_translate.return_value = 'HPHGDSSIYDAMVRMSQ'
        get_seq_diffs('CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA', geneToRef['PARC'])
        self.assertEqual(mock_six_frame_translate.call_args_list, [call('CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA', 1)])

    def test_update_GBS_Res_var(self):
        # ============== Test PARC ==================
        update_GBS_Res_var('PARC', ['Q17S'], GBS_Res_var)
        self.assertEqual(GBS_Res_var, {
            'PARC':'PARC-Q17S',
            'GYRA': '',
            '23S1': '',
            '23S3': '',
            'RPOBGBS-1': '',
            'RPOBGBS-2': '',
            'RPOBGBS-3': '',
            'RPOBGBS-4': ''
        })

        # ============== Test PARC with new variant ==================
        update_GBS_Res_var('PARC', ['Q18S'], GBS_Res_var)
        self.assertEqual(GBS_Res_var, {
            'PARC':'PARC-Q18S',
            'GYRA': '',
            '23S1': '',
            '23S3': '',
            'RPOBGBS-1': '',
            'RPOBGBS-2': '',
            'RPOBGBS-3': '',
            'RPOBGBS-4': ''
        })

        # ============== Test GYRA with no variant ==================
        update_GBS_Res_var('GYRA', [], GBS_Res_var)
        self.assertEqual(GBS_Res_var, {
            'PARC':'PARC-Q18S',
            'GYRA': 'GYRA',
            '23S1': '',
            '23S3': '',
            'RPOBGBS-1': '',
            'RPOBGBS-2': '',
            'RPOBGBS-3': '',
            'RPOBGBS-4': ''
        })

        # ============== Test 23S1 with two SNPs ==================
        update_GBS_Res_var('23S1', ['G1A', 'G34T'], GBS_Res_var)
        self.assertEqual(GBS_Res_var, {
            'PARC':'PARC-Q18S',
            'GYRA': 'GYRA',
            '23S1': '23S1-G1A,G34T',
            '23S3': '',
            'RPOBGBS-1': '',
            'RPOBGBS-2': '',
            'RPOBGBS-3': '',
            'RPOBGBS-4': ''
        })

    def test_update_drug_res_col_dict(self):
        # ============== Test PARC variant ==================
        update_drug_res_col_dict('PARC', ['Q17S'], drugRes_Col, geneToClass)
        self.assertEqual(drugRes_Col, {
            'TET': 'neg',
            'EC': 'neg',
            'FQ': 'PARC-Q17S',
            'OTHER': 'neg'
        })

        # ============== Test RPOBgbs-1 variant ==================
        update_drug_res_col_dict('RPOBGBS-1', ['F1G'], drugRes_Col, geneToClass)
        self.assertEqual(drugRes_Col, {
            'TET': 'neg',
            'EC': 'neg',
            'FQ': 'PARC-Q17S',
            'OTHER': 'RPOBGBS-1-F1G'
        })

    def test_get_seq_content(self):
        actual = get_seq_content(self.TEST_CONSENSUS_SEQ_FILE)
        self.assertEqual(actual, {
            '11__23S1__23S1-1__11': 'GTTACCCGCGACAGGACGGAAAGACCCCATGGAG',
            '12__23S3__23S3-3__12': 'CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG',
            '16__RPOBgbs__RPOBgbs-1__16': 'TTTGGTTCATCACAGCTGTCACAATTCATGGACCAACACAACCCTCTATCAGAATTGTCGCACAAACGCCGTCTCTCTGCCTTAGGACCTGGTGGTTTG',
            '17__RPOBgbs__RPOBgbs-2__17': 'GTTTCACAATTAGTCCGTTCTCCTGGTGTT',
            '18__RPOBgbs__RPOBgbs-3__18': 'TTTACAGTTGCACAAGCCAACTCTAAGCTTAACGAAGACGGTACATTTGCAGAAGAAATCGTTATGGGTCGTCATCAAGGTAATAACCAAGAGTTTCCTTCAAGCATT',
            '19__RPOBgbs__RPOBgbs-4__19': 'TTGATTGATCCAAAAGCACCATATGTTGGTACT',
            '5__GYRAGBS__GYRAGBS-1__5': 'GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATGGCACAATGGTGG',
            '7__PARCGBS__PARCGBS-1__7': 'CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA'
        })

    def test_get_gene_names_from_consensus(self):
        consensus_seq_dict = {
            '11__23S1__23S1-1__11': 'GTTACCCGCGACAGGACGGAAAGACCCCATGGAG',
            '12__23S3__23S3-3__12': 'CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG',
            '16__RPOBgbs__RPOBgbs-1__16': 'TTTGGTTCATCACAGCTGTCACAATTCATGGACCAACACAACCCTCTATCAGAATTGTCGCACAAACGCCGTCTCTCTGCCTTAGGACCTGGTGGTTTG',
            '17__RPOBgbs__RPOBgbs-2__17': 'GTTTCACAATTAGTCCGTTCTCCTGGTGTT',
            '18__RPOBgbs__RPOBgbs-3__18': 'TTTACAGTTGCACAAGCCAACTCTAAGCTTAACGAAGACGGTACATTTGCAGAAGAAATCGTTATGGGTCGTCATCAAGGTAATAACCAAGAGTTTCCTTCAAGCATT',
            '19__RPOBgbs__RPOBgbs-4__19': 'TTGATTGATCCAAAAGCACCATATGTTGGTACT',
            '5__GYRAGBS__GYRAGBS-1__5': 'GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATGGCACAATGGTGG',
            '7__PARCGBS__PARCGBS-1__7': 'CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA'
        }
        actual = get_gene_names_from_consensus(consensus_seq_dict)
        self.assertEqual(actual, ['PARC','GYRA','23S1','23S3','RPOBGBS-1','RPOBGBS-2','RPOBGBS-3','RPOBGBS-4'])

    @patch('bin.process_res_typer_results.get_seq_content')
    @patch('bin.process_res_typer_results.get_gene_names_from_consensus')
    @patch('bin.process_res_typer_results.get_seq_diffs')
    @patch('bin.process_res_typer_results.update_GBS_Res_var')
    @patch('bin.process_res_typer_results.update_drug_res_col_dict')
    def test_get_variants(self, mock_update_drug_res_col_dict, mock_update_GBS_Res_var, mock_get_seq_diffs, mock_get_gene_names_from_consensus, mock_get_seq_content):
        mock_get_seq_content.return_value = {
            '11__23S1__23S1-1__11': 'GTTACCCGCGACAGGACGGAAAGACCCCATGGAG',
            '12__23S3__23S3-3__12': 'CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG',
            '16__RPOBgbs__RPOBgbs-1__16': 'TTTGGTTCATCACAGCTGTCACAATTCATGGACCAACACAACCCTCTATCAGAATTGTCGCACAAACGCCGTCTCTCTGCCTTAGGACCTGGTGGTTTG',
            '17__RPOBgbs__RPOBgbs-2__17': 'GTTTCACAATTAGTCCGTTCTCCTGGTGTT',
            '18__RPOBgbs__RPOBgbs-3__18': 'TTTACAGTTGCACAAGCCAACTCTAAGCTTAACGAAGACGGTACATTTGCAGAAGAAATCGTTATGGGTCGTCATCAAGGTAATAACCAAGAGTTTCCTTCAAGCATT',
            '19__RPOBgbs__RPOBgbs-4__19': 'TTGATTGATCCAAAAGCACCATATGTTGGTACT',
            '5__GYRAGBS__GYRAGBS-1__5': 'GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATGGCACAATGGTGG',
            '7__PARCGBS__PARCGBS-1__7': 'CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA'
        }
        mock_get_gene_names_from_consensus.return_value = ['PARC','GYRA','23S1','23S3','RPOBGBS-1','RPOBGBS-2','RPOBGBS-3','RPOBGBS-4']
        mock_get_seq_diffs.return_value = ['Q17S']
        get_variants(self.TEST_CONSENSUS_SEQ_FILE)
        self.assertEqual(mock_get_seq_content.call_args_list, [call(self.TEST_CONSENSUS_SEQ_FILE)])
        self.assertEqual(mock_get_gene_names_from_consensus.call_args_list, [call(mock_get_seq_content.return_value)])
        self.assertEqual(mock_get_seq_diffs.call_args_list, [])
        self.assertEqual(mock_update_GBS_Res_var.call_args_list, [])
        self.assertEqual(mock_update_drug_res_col_dict.call_args_list, [])

    @patch('bin.process_res_typer_results.derive_presence_absence_targets')
    @patch('bin.process_res_typer_results.derive_presence_absence_targets_for_arg_res')
    @patch('lib.file_utils.FileUtils.create_output_contents')
    @patch('bin.process_res_typer_results.get_variants')
    @patch('lib.file_utils.FileUtils.write_output')
    def test_run(self, mock_write_output, mock_get_variants, mock_create_output_contents, mock_derive_presence_absence_targets_for_arg_res, mock_derive_presence_absence_targets):
        args = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_gbs_fullgenes', '--srst2_gbs_consensus', 'srst2_gbs_consensus',
            '--srst2_other_fullgenes', 'srst2_argannot_fullgenes', 'srst2_resfinder_fullgenes',
            '--min_read_depth', '30', '--output_prefix', 'output'])
        mock_create_output_contents.return_value = 'foobar'
        run(args)
        self.assertEqual(mock_derive_presence_absence_targets.call_args_list, [call(args.srst2_gbs_fg_output)])
        self.assertEqual(mock_derive_presence_absence_targets_for_arg_res.call_args_list, [call(args.srst2_other_fg_output)])
        mock_create_output_contents.assert_has_calls([
            call(Res_Targets),
            call(GBS_Res_var),
            call(drugRes_Col)
        ], any_order = False)
        self.assertEqual(mock_get_variants.call_args_list, [call(args.srst2_gbs_cs_output)])
        mock_write_output.assert_has_calls([
            call(ANY, args.output + '_res_incidence.txt'),
            call(ANY, args.output + "_res_gbs_variants.txt"),
            call(ANY, args.output + "_res_alleles_variants.txt")
        ], any_order = False)

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_gbs_fullgenes', '--srst2_gbs_consensus', 'srst2_gbs_consensus',
            '--srst2_other_fullgenes', 'srst2_argannot_fullgenes', 'srst2_resfinder_fullgenes',
            '--min_read_depth', '30.0', '--output_prefix', 'output'])
        self.assertEqual(actual,
                         argparse.Namespace(srst2_gbs_fg_output='srst2_gbs_fullgenes',
                                            srst2_gbs_cs_output='srst2_gbs_consensus',
                                            srst2_other_fg_output=['srst2_argannot_fullgenes','srst2_resfinder_fullgenes'],
                                            min_depth = 30.0,
                                            output='output'))

    @patch('bin.process_res_typer_results.get_arguments')
    @patch('bin.process_res_typer_results.run')
    def test_main(self, mock_run, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        main()
        self.assertEqual(ANY, [call(args)])
