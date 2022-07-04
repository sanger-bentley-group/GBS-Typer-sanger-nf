#!/usr/bin/env python3
import argparse
import sys
import os
import re
import glob
import subprocess
import json
import pandas as pd
from collections import defaultdict
from lib.six_frame_translation import six_frame_translate, extract_frame_aa, codon2aa
from lib.file_io import get_seq_content
from lib.file_utils import FileUtils
from lib.get_headers import read_header_json

class nSeq(str): # Nucleotide sequence
    pass


class aSeq(str): # Amino acid sequence
    pass


# Drug Class Resistance dictionary
drugRes_Col = {
    'AG': 'neg',
    'TET': 'neg',
    'EC': 'neg',
    'FQ': 'neg',
    'OTHER': 'neg',
}

# Gene to Drug Class Resistance lookup dictionary
geneToClass = {
    'AAC6APH2': 'AG',
    'AADECC': 'AG',
    'ANT6IA3KF864551': 'AG',
    'ANT6IA': 'AG',
    'APH3III': 'AG',
    'APH3OTHER': 'AG',
    'CATPC194': 'OTHER',
    'CATQ': 'OTHER',
    'ERMA': 'EC',
    'ERMB': 'EC',
    'ERMT': 'EC',
    'FOSA': 'OTHER',
    'LNUB': 'EC',
    'LNUC': 'EC',
    'LSAC': 'EC',
    'LSAE': 'EC',
    'MEFA': 'EC',
    'MPHC': 'EC',
    'MSRA': 'OTHER',
    'MSRD': 'OTHER',
    'TETB': 'TET',
    'TETL': 'TET',
    'TETM': 'TET',
    'TETO32O': 'TET',
    'TETOW32OWO': 'TET',
    'TETOW32O': 'TET',
    'TETOWO': 'TET',
    'TETOW': 'TET',
    'TETO': 'TET',
    'TETSM': 'TET',
    'TETS': 'TET',
    'TETW32O': 'TET',
    'TETW4FN396364': 'TET',
    'SUL2': 'OTHER',
    'PARC': 'FQ',
    'GYRA': 'FQ',
    '23S1': 'EC',
    '23S3': 'EC',
    'RPOBGBS-1': 'OTHER',
    'RPOBGBS-2': 'OTHER',
    'RPOBGBS-3': 'OTHER',
    'RPOBGBS-4': 'OTHER',
}

# Other Resistance Targets dictionary
Res_Targets = {
    'AAC6APH2': 'neg',
    'AADECC': 'neg',
    'ANT6IA3KF864551': 'neg',
    'ANT6IA': 'neg',
    'APH3III': 'neg',
    'APH3OTHER': 'neg',
    'CATPC194': 'neg',
    'CATQ': 'neg',
    'ERMA': 'neg',
    'ERMB': 'neg',
    'ERMT': 'neg',
    'FOSA': 'neg',
    'LNUB': 'neg',
    'LNUC': 'neg',
    'LSAC': 'neg',
    'LSAE': 'neg',
    'MEFA': 'neg',
    'MPHC': 'neg',
    'MSRA': 'neg',
    'MSRD': 'neg',
    'SUL2': 'neg',
    'TETB': 'neg',
    'TETL': 'neg',
    'TETM': 'neg',
    'TETO32O': 'neg',
    'TETOW32OWO': 'neg',
    'TETOW32O': 'neg',
    'TETOWO': 'neg',
    'TETOW': 'neg',
    'TETO': 'neg',
    'TETSM': 'neg',
    'TETS': 'neg',
    'TETW32O': 'neg',
    'TETW4FN396364': 'neg'
}

# GBS Resistance Targets dictionary
GBS_Res_Targets = {
    'GYRA': 'neg',
    'PARC': 'neg',
    '23S1': 'neg',
    '23S3': 'neg',
    'RPOBGBS-1': 'neg',
    'RPOBGBS-2': 'neg',
    'RPOBGBS-3': 'neg',
    'RPOBGBS-4': 'neg',
}

# GBS Gene Resistance Variants dictionary
GBS_Res_var = {
    'GYRA_SNP': '', #10
    'PARC_SNP': '', #14
    '23S1_SNP': '', #0
    '23S3_SNP': '', #1
    'RPOBGBS-1_SNP': '', #6
    'RPOBGBS-2_SNP': '', #7
    'RPOBGBS-3_SNP': '', #8
    'RPOBGBS-4_SNP': '', #9,
}

# Reference sequence dictionary
geneToRef = defaultdict(lambda: '')
geneToRef.update({
    'PARC': aSeq('HPHGDSSIYDAMVRMSQ'),
    'GYRA': aSeq('VMGKYHPHGDSSIYEAMVRMAQWW'),
    '23S1': nSeq('GTTACCCGCGACAGGACGGAAAGACCCCATGGAG'),
    '23S3': nSeq('CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG'),
    'RPOBGBS-1': aSeq('FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL'),
    'RPOBGBS-2': aSeq('VSQLVRSPGV'),
    'RPOBGBS-3': aSeq('FTVAQANSKLNEDGTFAEEIVMGRHQGNNQEFPSSI'),
    'RPOBGBS-4': aSeq('LIDPKAPYVGT')
})

# Gene to target name lookup dictionary
geneToTargetSeq = defaultdict(lambda: '')
geneToTargetSeq.update({
    'PARC': '7__PARCGBS__PARCGBS-1__7',
    'GYRA': '5__GYRAGBS__GYRAGBS-1__5',
    '23S1': '11__23S1__23S1-1__11',
    '23S3': '12__23S3__23S3-3__12',
    'RPOBGBS-1': '16__RPOBgbs__RPOBgbs-1__16',
    'RPOBGBS-2': '17__RPOBgbs__RPOBgbs-2__17',
    'RPOBGBS-3': '18__RPOBgbs__RPOBgbs-3__18',
    'RPOBGBS-4': '19__RPOBgbs__RPOBgbs-4__19'
})

snpOffset = defaultdict(lambda: 0)
snpOffset.update({
    "GYRA": 70,
    "PARC": 73
})

geneAlleleDict = defaultdict(lambda: [])

EOL_SEP = "\n"


def update_presence_absence_target(gene, allele, depth, gbs_res_target_dict):
    """Update presence/absence for GBS Targets dictionary"""
    if depth >= MIN_DEPTH:

        for gene_name in gbs_res_target_dict.keys():
            if re.search(gene_name, allele):

                gbs_res_target_dict[gene_name] = "pos"


def derive_presence_absence_targets(input_file, GBS_Res_Targets):
    """Find gene presence/absence for the GBS resistance database"""
    try:
        with open(input_file, 'r') as fd:
            # Skip header row
            next(fd)

            # Process file lines
            for line in fd:
                fields = line.split('\t')
                gene = fields[2]
                allele = fields[3]
                depth = float(fields[5])
                update_presence_absence_target(gene, allele, depth, GBS_Res_Targets)
    except IOError:
        print('Cannot open {}.'.format(input_file))


def update_presence_absence_target_for_arg_res(gene, allele, depth, drug_res_col_dict, res_target_dict, gene_allele_dict):
    """Update presence/absence for Other Resistance Targets dictionary"""
    if depth >= MIN_DEPTH:

        other = 1
        for gene_name in res_target_dict.keys():
            if re.search(gene_name, "".join(re.split("[^a-zA-Z0-9]*", allele)).upper()):
                other = 0

                if allele not in gene_allele_dict.keys():
                    gene_allele_dict[allele] = gene_name

                if res_target_dict[gene_name] == "neg":
                    res_target_dict[gene_name] = "pos"

                drugCat = geneToClass[gene_name]
                if drug_res_col_dict[drugCat] == "neg":
                    drug_res_col_dict[drugCat] = gene + '[' + allele + ']'
                else:
                    drug_res_col_dict[drugCat] = drug_res_col_dict[drugCat] + ':' + gene + '[' + allele + ']'

                break

        if other:
            if drug_res_col_dict["OTHER"] == "neg":
                drug_res_col_dict["OTHER"] = gene + '[' + allele + ']'
            else:
                drug_res_col_dict["OTHER"] = drug_res_col_dict["OTHER"] + ":" + gene + '[' + allele + ']'


def clear_arg_res(res_target_dict):
    for key in res_target_dict:
        res_target_dict[key] = ''


def derive_presence_absence_targets_for_arg_res(input_files, drugRes_Col, Res_Targets):
    """Find gene presence/absence for other resistance databases"""
    for input_file in input_files:
        if os.stat(input_file).st_size != 0:
            try:
                with open(input_file, 'r') as fd:
                    # Skip header row
                    next(fd)
                    # Process file lines
                    for line in fd:
                        fields = line.split('\t')
                        gene = fields[2]
                        allele = fields[3]
                        depth = float(fields[5])
                        update_presence_absence_target_for_arg_res(gene, allele, depth, drugRes_Col, Res_Targets, geneAlleleDict)
            except IOError:
                print('Cannot open {}.'.format(input_file))
        else:
            clear_arg_res(Res_Targets)


def find_mismatches(seq_diffs, query_Seq, ref_Seq, snp_Offset):
    """Find mismatches between query and reference sequences"""
    for resi in range(len(query_Seq)):
        if query_Seq[resi] != ref_Seq[resi]:
            seq_diffs.append(ref_Seq[resi] + str(resi+1+snp_Offset) + query_Seq[resi])
    return seq_diffs


def get_seq_diffs(query_Seq, ref_Seq, snp_Offset):
    """Get SNP variants"""
    if type(ref_Seq) == aSeq:
        query_Seq = six_frame_translate(query_Seq, 1)
    seq_diffs = []
    if query_Seq != ref_Seq:
        seq_diffs = find_mismatches(seq_diffs, query_Seq, ref_Seq, snp_Offset)
    return seq_diffs


def update_GBS_Res_var(gene_name, seq_diffs, bin_res_arr):
    """Update GBS Gene Resistance Variants dictionary with gene variants"""
    if seq_diffs:
        bin_res_arr[gene_name + '_SNP'] = ','.join(seq_diffs)
    else:
        bin_res_arr[gene_name + '_SNP'] = '*'


def update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, geneToClass):
    """Update Drug Resistance Class dictionary with GBS gene variants"""
    drugClass = geneToClass[gene_name]
    if seq_diffs:
        gene_var = gene_name + '-' + ','.join(seq_diffs)
    else:
        gene_var = gene_name
    if drugRes_Col[drugClass] == 'neg':
        drugRes_Col[drugClass] = gene_var
    else:
        new_value = drugRes_Col[drugClass] + ':' + gene_var
        drugRes_Col[drugClass] = new_value


def get_gene_names_from_consensus(consensus_seq_dict):
    """Get the gene name identifiers from the consensus sequence IDs"""
    gene_names = []
    for gene_name in geneToTargetSeq.keys():
        if geneToTargetSeq[gene_name] in consensus_seq_dict.keys():
            gene_names.append(gene_name)
    return gene_names


def get_variants(consensus_seqs):
    """Get resistance gene variants from freebayes consensus GBS sequences"""
    consensus_seq_dict = get_seq_content(consensus_seqs)
    gene_names = get_gene_names_from_consensus(consensus_seq_dict)
    for gene_name in gene_names:
        if GBS_Res_Targets[gene_name] == "pos" and geneToTargetSeq[gene_name] and geneToRef[gene_name]:
            seq_diffs = get_seq_diffs(consensus_seq_dict[geneToTargetSeq[gene_name]], geneToRef[gene_name], snpOffset[gene_name])
            update_GBS_Res_var(gene_name, seq_diffs, GBS_Res_var)
            update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, geneToClass)


def run(args):

    # Set minimum read depth
    global MIN_DEPTH
    MIN_DEPTH = args.min_depth

    # Get presence/absence of genes
    derive_presence_absence_targets(args.srst2_gbs_fg_output, GBS_Res_Targets)

    if args.srst2_other_fg_output is not None:
        derive_presence_absence_targets_for_arg_res(args.srst2_other_fg_output, drugRes_Col, Res_Targets)
        GBS_Res_Targets.update(Res_Targets)

    inc_out = FileUtils.create_output_contents(GBS_Res_Targets)

    # Get variants
    get_variants(args.srst2_gbs_cs_output)
    var_out = FileUtils.create_output_contents(GBS_Res_var)

    # Get alleles for all drug classes
    allele_out = FileUtils.create_output_contents(drugRes_Col)

    # Write allele accessions/gene contents
    header_dict = read_header_json(args.headers)
    content = ""
    id = args.output.split('/')[len(args.output.split('/'))-1]
    for key, value in geneAlleleDict.items():
        if value in header_dict["combine_all"]:
            print(value)
            row_name = header_dict["combine_all"][value]
            content = f'{content}{id}\t{row_name}\t{key}\n'

    try:
        output_filename = args.output + "_res_alleles_accessions.txt"
        with open(output_filename, 'w') as out:
            out.write(content)
    except IOError:
        print('Cannot open filename starting "{}"'.format(output_filename))
        raise

    # Write incidence output
    FileUtils.write_output(inc_out, args.output + '_res_incidence.txt')
    # Write gbs variant output
    FileUtils.write_output(var_out, args.output + "_res_gbs_variants.txt")
    # Write allele output
    FileUtils.write_output(allele_out, args.output + "_res_alleles_variants.txt")


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify SRST2 sequence typing output files.')
    parser.add_argument('--srst2_gbs_fullgenes', dest='srst2_gbs_fg_output', required=True,
                        help='Input SRST2 fullgenes output for the GBS reference database.')
    parser.add_argument('--srst2_gbs_consensus', dest='srst2_gbs_cs_output', required=True,
                        help='Input freebayes consensus sequence output for the GBS reference database.')
    parser.add_argument('--srst2_other_fullgenes', dest='srst2_other_fg_output', required=False,
                        help='Input SRST2 fullgenes outputs for other references databases.',
                        nargs = '*')
    parser.add_argument('--min_read_depth', dest='min_depth', required=True, type=float, default=30,
                        help = 'Minimum read depth where mappings with fewer reads are excluded. Default: 30.')
    parser.add_argument('--headers', dest='headers', required=True,
                        help='JSON file of expected headers.')
    parser.add_argument('--output_prefix', dest='output', required=True,
                        help='Output prefix of filename.')

    return parser


def main():
    args = get_arguments().parse_args()
    run(args)


if __name__ == "__main__":
    sys.exit(main())
