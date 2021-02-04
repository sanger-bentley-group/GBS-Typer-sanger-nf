from collections import defaultdict


def get_seq_content(seq_file):
    """Return sequence file content in dictionary"""
    seq_dict = defaultdict(lambda: '')
    try:
        with open(seq_file, 'r') as fd:
            for line in fd:
                if line[0] == '>':
                    seq_name = line.split('>')[1].split('\n')[0]
                else:
                    if seq_dict[seq_name] == '':
                        seq_dict[seq_name] = line.split('\n')[0]
                    else:
                        tmp_seq = seq_dict[seq_name] + line.split('\n')[0]
                        seq_dict[seq_name] = tmp_seq
    except IOError:
        print('Cannot open {}.'.format(seq_file))

    return seq_dict


def write_seq_dict(seq_dict, output_filename):
    try:
        with open(output_filename, 'w') as out:
            for header in seq_dict:
                out.write('>' + header + '\n' + seq_dict[header] + '\n')

    except IOError:
        print('Cannot open {}.'.format(output_filename))
