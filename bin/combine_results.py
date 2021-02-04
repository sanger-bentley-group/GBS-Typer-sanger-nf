#!/usr/bin/env python3
from __future__ import print_function
import argparse
import sys
from lib.file_utils import FileUtils


def get_content_with_id(id, file):
    """ Include ID string with content of file """
    content = 'ID' + '\t'
    count = 0
    with open(file, 'r') as f:
        for line in f:
            if count:
                content += id + '\t' + line
            else:
                content += line
            count += 1
    return content


def get_sero_res_contents(id, sero_file, res_file):
    """ Merge serotyping and resistance typing content"""
    header = 'ID' + '\t' + 'Serotype' + '\t'
    type = ''
    res_incidence = ''
    with open(sero_file, 'r') as sero:
        next(sero) # Skip header row
        count = 0
        for line in sero:
            if count:
                type = type + ";" + line.split('\t')[2]
            else:
                type = id + '\t' + line.split('\t')[2]
            count += 1

    with open(res_file, 'r') as res:
        count = 0
        for line in res:
            if count:
                results = line.split('\n')[0].split('\t')
                results_out = []
                for result in results:
                    if result == 'pos':
                        results_out.append('+')
                    elif result == 'neg':
                        results_out.append('-')
                res_incidence = '\t'.join(results_out)
            else:
                header += line
            count += 1

    return(header + type + '\t' + res_incidence + '\n')


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
    subparser_pbp_typing.add_argument('--pbp_existing_allele_results', '-p', dest='pbp_allele', required=True,
                        help='Input surface typing incidence results file.')
    subparser_pbp_typing.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    subparser_pbp_typing.set_defaults(which='pbp_typer')

    return parser


def main():
    parser = get_arguments()
    args = parser.parse_args()

    if args.which == "sero_res":
        # Merge serotyping and resistance typing results (including ID)
        sero_res_output_lines = get_sero_res_contents(args.id, args.sero, args.inc)
        FileUtils.write_output(sero_res_output_lines, args.output + "_sero_res_incidence.txt")

        # Add ID to alleles from resistance typing results
        res_alleles_output_lines = get_content_with_id(args.id, args.alleles)
        FileUtils.write_output(res_alleles_output_lines, args.output + "_id_alleles_variants.txt")

        # Add ID to variants from resistance typing results
        res_variants_output_lines = get_content_with_id(args.id, args.variants)
        FileUtils.write_output(res_variants_output_lines, args.output + "_id_variants.txt")

    elif args.which == "surface_typer":
        # Add ID to surface typing incidence results
        if args.surface_inc:
            surface_protein_incidence_output_lines = get_content_with_id(args.id, args.surface_inc)
            FileUtils.write_output(surface_protein_incidence_output_lines, args.output + "_surface_protein_incidence.txt")

        # Add ID to surface typing variants results
        if args.surface_variants:
            surface_protein_variants_output_lines = get_content_with_id(args.id, args.surface_variants)
            FileUtils.write_output(surface_protein_variants_output_lines, args.output + "_surface_protein_variants.txt")

    elif args.which == "pbp_typer":
        # Add ID to PBP typer existing allele results
        if args.pbp_allele:
            pbp_typer_output_lines = get_content_with_id(args.id, args.pbp_allele)
            FileUtils.write_output(pbp_typer_output_lines, args.output + "_existing_PBP_allele.txt")

    else:
        print("ERROR: Please specify a valid option.")
        parser.print_help()


if __name__ == "__main__":
    sys.exit(main())
