#!/usr/bin/env python3
from collections import defaultdict
import argparse, sys
from lib.seq_data import SeqData, BlastData


class FragmentData():
    def __init__(self):
        self._fragment_positions = {}

    def get_start_end_positions(self, best_hits, seq_lengths, alignment_length_threshold, identity_threshold):
        for allele, stat in best_hits.items():
            if float(best_hits[allele][1]) >= identity_threshold * 100 and float(best_hits[allele][2]) >= alignment_length_threshold * seq_lengths[allele]:
                self.calculate_start_end_positions(allele, best_hits[allele], seq_lengths[allele])

    def calculate_start_end_positions(self, allele, stats, seq_length):
        if int(stats[8]) > int(stats[7]):
            frag_start = int(stats[7]) - int(stats[5])
            frag_end = (seq_length - int(stats[6])) + int(stats[8])
            self._fragment_positions[allele] = (stats[0], str(frag_start), str(frag_end), 'forward', '1', '+')
        elif int(stats[8]) < int(stats[7]):
            frag_start = int(stats[7]) + int(stats[5]) - 1
            frag_end = int(stats[8]) - (seq_length - int(stats[6])) - 1
            self._fragment_positions[allele] = (stats[0], str(frag_end), str(frag_start), 'reverse', '1', '-')

    def get_data(self):
        return self._fragment_positions

    def write_start_end_positions(self, output_prefix):
        if self._fragment_positions.items():
            for allele, frag_stats in self._fragment_positions.items():
                output_filename = output_prefix + allele + '.bed'
                try:
                    with open(output_filename, 'w') as out:
                        out.write('\t'.join(frag_stats) + '\n')
                except IOError:
                    print('Cannot open {}.'.format(output_filename))


def get_arguments():
    parser = argparse.ArgumentParser(description='Output the start and end positions of b-lactam genes from contigs in a BED file.')
    parser.add_argument('--blast_out_file', '-b', dest='blast_out', required=True,
                        help='Input BLAST results file.', type = str)
    parser.add_argument('--query_fasta', '-f', dest='fasta_qu', required=True,
                        help='Fasta file query.', type = str)
    parser.add_argument('--frac_align_len_threshold', '-fa', dest='frac_align', required=False,
                        help='Fraction of alignment length threshold.', type = float, default=0.5)
    parser.add_argument('--frac_identity_threshold', '-fi', dest='frac_ident', required=False,
                        help='Fraction of identity threshold.', type = float, default=0.5)
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output BED file.', type = str)
    return parser


def check_arguments(args):
    if not (0 <= args.frac_align <= 1):
        raise ValueError("Invalid frac_align_len_threshold value. Value must be between 0 and 1.")

    if not (0 <= args.frac_ident <= 1):
        raise ValueError("Invalid frac_identity_threshold value. Value must be between 0 and 1.")


def main():
    args = get_arguments().parse_args()

    check_arguments(args)

    # Load blast data
    blast_data = BlastData(args.blast_out)

    # Get blast best hit
    best_blast_hits = blast_data.get_best_hit()

    # Load sequence data
    seq_data = SeqData(args.fasta_qu)

    # Get lengths of reference sequences
    seq_lengths = seq_data.calculate_seq_length()

    # Get start and end positions of fragments
    fragment_data = FragmentData()
    fragment_data.get_start_end_positions(best_blast_hits, seq_lengths, args.frac_align, args.frac_ident)

    # Write fragments
    fragment_data.write_start_end_positions(args.output)


if __name__ == "__main__":
    sys.exit(main())
