#!/usr/bin/env python3
from __future__ import print_function
import argparse
import sys
from lib.file_utils import FileUtils
from lib.get_headers import read_header_json
import pandas as pd
import json


def get_content(file):
    try:
        df = pd.read_csv(file, sep="\t")
    except:
        df = pd.DataFrame()

    return df


def create_model_df(headers, id_df):
    model_df = pd.DataFrame(columns=headers, index = [0])
    model_df = id_df.merge(model_df, how="inner", left_index=True, right_index=True)

    return model_df


def merge_dfs(df1, df2):
    df1 = df1.combine_first(df2)

    return df1


def create_df(headers: list, id_df: pd.DataFrame, files: list):
    output_df = create_model_df(headers, id_df)

    for file in files:
        content_df = get_content(file)
        output_df = merge_dfs(output_df, content_df)

    headers = id_df.columns.to_list() + headers

    return output_df.loc[:,headers]


def rename_columns(df, header_dict: dict, id_df: pd.DataFrame):
    df.rename(columns = header_dict, inplace = True)
    headers = id_df.columns.to_list() + list(header_dict.values())
    df = df.loc[:,headers]

    return df


def get_arguments():
    """Parse allowed argument combinations"""
    parser = argparse.ArgumentParser(description='Combine sample results for a specified pipeline.')
    subparsers = parser.add_subparsers(title='Available commands', help='', metavar='')

    subparser_sero_res = subparsers.add_parser(
        'sero_res',
        help='',
        description='Combine resistance gene typer and serotyper results.',
    )

    subparser_sero_res.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID.')
    subparser_sero_res.add_argument('--headers', '-t', dest='headers', required=True,
                        help='JSON file of expected headers.')
    subparser_sero_res.add_argument('--serotyper_results', '-s', dest='sero', required=True,
                        help='Input SeroType results tab file.')
    subparser_sero_res.add_argument('--res_incidence_results', '-r', dest='inc', required=True,
                        help='Input resistance typing incidence results file.')
    subparser_sero_res.add_argument('--res_alleles_results', '-a', dest='alleles', required=True,
                        help='Input resistance typing alleles results file.')
    subparser_sero_res.add_argument('--res_variants_results', '-v', dest='variants', required=True,
                        help='Input resistance typing variants results file.')
    subparser_sero_res.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    subparser_sero_res.set_defaults(which='sero_res')

    subparser_surface_typing = subparsers.add_parser(
        'surface_typer',
        help='',
        description='Combine surface protein typing results.',
    )

    subparser_surface_typing.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID.')
    subparser_surface_typing.add_argument('--headers', '-t', dest='headers', required=True,
                        help='JSON file of expected headers.')
    subparser_surface_typing.add_argument('--surface_incidence_results', '-x', dest='surface_inc', required=True,
                        help='Input surface typing incidence results file.')
    subparser_surface_typing.add_argument('--surface_variants_results', '-y', dest='surface_variants', required=True,
                        help='Input surface typing variants results file.')
    subparser_surface_typing.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    subparser_surface_typing.set_defaults(which='surface_typer')

    subparser_pbp_typing = subparsers.add_parser(
        'pbp_typer',
        help='',
        description='Combine PBP typer results.',
    )

    subparser_pbp_typing.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID.')
    subparser_pbp_typing.add_argument('--headers', '-t', dest='headers', required=True,
                        help='JSON file of expected headers.')
    subparser_pbp_typing.add_argument('--pbp_existing_allele_results', '-p', dest='pbp_allele', required=True,
                        help='Input surface typing incidence results file.')
    subparser_pbp_typing.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    subparser_pbp_typing.set_defaults(which='pbp_typer')

    subparser_combine_all = subparsers.add_parser(
        'combine_all',
        help='',
        description='Combine all results.'
    )

    subparser_combine_all.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID.')
    subparser_combine_all.add_argument('--headers', '-t', dest='headers', required=True,
                        help='JSON file of expected headers.')
    subparser_combine_all.add_argument('--serotyper_results', '-s', dest='sero', required=True,
                        help='Input SeroType results tab file.')
    subparser_combine_all.add_argument('--res_incidence_results', '-r', dest='inc', required=True,
                        help='Input resistance typing incidence results file.')
    subparser_combine_all.add_argument('--res_variants_results', '-v', dest='variants', required=True,
                        help='Input resistance typing variants results file.')
    subparser_combine_all.add_argument('--mlst_allelic_frequency_results', '-m', dest='mlst', required=True,
                        help='Input SRST2 results file of MLST allelic frequency.')
    subparser_combine_all.add_argument('--surface_incidence_results', '-x', dest='surface_inc', required=True,
                        help='Input surface typing incidence results file.')
    subparser_combine_all.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    subparser_combine_all.set_defaults(which='combine_all')

    return parser


def main():
    parser = get_arguments()
    args = parser.parse_args()

    header_dict = read_header_json(args.headers)
    id_df = pd.DataFrame(args.id, columns=header_dict["id"], index = [0])

    if args.which == "sero_res":
        # Merge serotyping and resistance typing results (including ID)
        df_sero_res = create_df(header_dict["sero_res"], id_df, [args.sero, args.inc])
        FileUtils.write_pandas_output(df_sero_res, args.output + "_sero_res_incidence.txt")

        # Add ID to alleles from resistance typing results
        df_res_alleles = create_df(header_dict["res_alleles"], id_df, [args.alleles])
        FileUtils.write_pandas_output(df_res_alleles, args.output + "_id_alleles_variants.txt")

        # Add ID to variants from resistance typing results
        df_gbs_res_variants = create_df(header_dict["gbs_res_variants"], id_df, [args.variants])
        FileUtils.write_pandas_output(df_gbs_res_variants, args.output + "_id_variants.txt")

    elif args.which == "surface_typer":
        # Add ID to surface typing incidence results
        if args.surface_inc:
            df_surface_inc = create_df(header_dict["surface_inc"], id_df, [args.surface_inc])
            FileUtils.write_pandas_output(df_surface_inc, args.output + "_surface_protein_incidence.txt")

        # Add ID to surface typing variants results
        if args.surface_variants:
            df_surface_variants = create_df(header_dict["surface_variants"], id_df, [args.surface_variants])
            FileUtils.write_pandas_output(df_surface_variants, args.output + "_surface_protein_variants.txt")

    elif args.which == "pbp_typer":
        # Add ID to PBP typer existing allele results
        if args.pbp_allele:
            df_pbp_allele = create_df(header_dict["pbp_allele"], id_df, [args.pbp_allele])
            FileUtils.write_pandas_output(df_pbp_allele, args.output + "_existing_PBP_allele.txt")

    elif args.which == "combine_all":
        # Combine ID, serotyping and resistance typing incidence, resistance typing variants, MLST type and allelic frequency, surface protein incidence
        df_combine_all = create_df(list(header_dict["combine_all"].keys()), id_df, [args.sero, args.inc, args.variants, args.mlst, args.surface_inc])
        df_combine_all = rename_columns(df_combine_all, header_dict["combine_all"], id_df)
        FileUtils.write_pandas_output(df_combine_all, args.output + '_id_combined_output.txt')

    else:
        print("ERROR: Please specify a valid option.")
        parser.print_help()


if __name__ == "__main__":
    sys.exit(main())
