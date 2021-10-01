#!/usr/bin/env python3
import argparse
import sys

def make_gene_dict(input_file, depth_threshold):
    """Get features from SRST2 input file into dictionary depending on read depth threshold"""
    gene_dict = dict()
    with open(input_file, 'r') as fg_file:
        next(fg_file) # Skip header row
        for line in fg_file:
            feature = line.split('\t')
            if float(feature[5]) > depth_threshold:
                gene_dict[feature[2]] = feature[2:-1]
    return gene_dict


def write_outfile(gene_dict, out_file):
    """Write serotype, match type status and average read depth to output file"""
    matched_alleles = []
    match_type = []
    serotype = []
    avgdepth = []
    for key, values in gene_dict.items():
        status = 'identical'
        if values[4] != '':
            status = 'imperfect'
        matched_alleles = matched_alleles + [values[1]]
        match_type = match_type + [values[0] + '=' + status]
        serotype = serotype + [values[0]]
        avgdepth = avgdepth + [values[3]]
    with open(out_file, 'w') as out:
        out.write('Matched_Allele'+'\t'+'Match_Type'+'\t'+'Serotype'+'\t'+'AvgDepth'+'\n'+'/'.join(matched_alleles)+'\t'+'/'.join(match_type)+'\t'+'/'.join(serotype)+'\t'+'/'.join(avgdepth)+'\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify fullgenes output of SRST2.')
    parser.add_argument('--srst2_output', '-s', dest='id', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--sero_db', '-b', dest='db', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output filename.')
    parser.add_argument('--min_read_depth', '-d', dest='depth', default = 30, type=float,
                        help='Minimum read depth where mappings with fewer reads are excluded. Default: 30.')

    return parser


def main():
    args = get_arguments().parse_args()
    db_name = ' '.join(args.db.split('.')[:-1])

    # Specift fullgenes file path from ID and database name
    fullgenes_file = args.id + '__fullgenes__' + db_name + '__results.txt'

    # Get feature dictionary
    gene_dict = make_gene_dict(fullgenes_file, args.depth)

    # Write tab-delimited output file with serotype features
    write_outfile(gene_dict, args.output)


if __name__ == "__main__":
    sys.exit(main())
