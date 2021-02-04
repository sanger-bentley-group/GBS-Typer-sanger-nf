#!/usr/bin/env python3
import argparse, sys
from lib.seq_data import SeqData
from lib.file_io import write_seq_dict

def get_arguments():
    parser = argparse.ArgumentParser(description='Output the start and end positions of b-lactam genes from contigs in a BED file.')
    parser.add_argument('--blactam_fasta', '-f', dest='fasta', required=True,
                        help='Input FASTA file of beta-lactam genes.', type = str)
    parser.add_argument('--output_file', '-o', dest='output', required=True,
                        help='Output file prefix of translated beta-lactam sequences.', type = str)
    return parser


def main():
    args = get_arguments().parse_args()

    # Get sequence data for PBP genes
    pbp_genes = SeqData(args.fasta)

    # Translate and write multiple FAA files of translated amino acid sequences
    pbp_genes.translate_content(1)

    # Write sequence file
    write_seq_dict(pbp_genes.get_data(), args.output)


if __name__ == "__main__":
    sys.exit(main())
