#!/usr/bin/env python
import argparse, sys, re

def get_targets(targets_file):
    targets = []
    with open(targets_file, 'r') as txt:
        for line in txt:
            targets.append(line.split('\n')[0])
    return targets


def in_line(line, target):
    if (bool(re.search('^@HD|^\@SQ.*' + target + '|^@PG', line))==True):
        return True
    elif (line.split('\t')[2] == target):
        return True
    else:
        return False


def write_sam_file(sam_file, target, id):
    with open('CHECK_' + target + '_' + id + '_seq.sam', 'w') as out:
        with open(sam_file, 'r') as sam:
            for line in sam:
                if in_line(line, target):
                    out.write(line)


def write_target_sam_files(targets, sam_file, id):
    for target in targets:
        write_sam_file(sam_file, target, id)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get targets from sam file.')
    parser.add_argument('--sam_file', '-s', dest='sam', required=True,
                        help='Input sam file.')
    parser.add_argument('--target_file', '-t', dest='target', required=True,
                        help='Input target text file.')
    parser.add_argument('--id', '-i', dest='id', required=True,
                        help='Read ID.')
    return parser


def main():
    args = get_arguments().parse_args()
    targets = get_targets(args.target)
    write_target_sam_files(targets, args.sam, args.id)


if __name__ == "__main__":
    sys.exit(main())