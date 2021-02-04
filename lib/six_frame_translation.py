#!/usr/bin/env python3
import re
from Bio.Seq import Seq

def codon2aa(codon):
    """Translate codons to amino acids"""

    codon = codon.upper()

    codon_dict = {
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
    }

    try:
        result = codon_dict[codon]
    except KeyError:
        if re.search(r"GC.", codon):
            result = 'A'
        elif re.search(r"GG.", codon):
            result = 'G'
        elif re.search(r"CC.", codon):
            result = 'P'
        elif re.search(r"AC.", codon):
            result = 'T'
        elif re.search(r"GT.", codon):
            result = 'V'
        elif re.search(r"CG.", codon):
            result = 'R'
        elif re.search(r"TC.", codon):
            result = 'S'
        else:
            result = 'x'
            print("Bad codon " + codon + "!!")

    return result


def extract_frame_aa(sequence, frame):
    """
    :param sequence: dna bases
    :param frame: frame number to extract
    :return: Amino acid translates frame
    """

    bases = sequence
    if frame > 3:
        seq = Seq(sequence)
        bases = str(seq.reverse_complement())
        frame -= 3

    i = frame - 1
    protein = ""
    while i < len(bases)-2:
        protein += codon2aa(bases[i:i+3])
        i += 3

    return protein


def six_frame_translate(seq_input, frame):
    """
    :param seq_input: Fasta feature including id and sequence lines
    :param frame: Codon number:
    :return: protein translation
    """

    if frame < 1 or frame > 6:
        raise IndexError("Frame number argument is out of bounds: " + str(frame))

    lines = seq_input.splitlines()
    dna = ""
    for line in lines:
        if line.startswith('>'):
            continue
        else:
            dna += line.strip()

    return extract_frame_aa(dna, frame)
