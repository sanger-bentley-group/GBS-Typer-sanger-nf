#!/usr/bin/env python3
"""Common classes containing and processing FASTA and BLAST output data"""
from collections import defaultdict
from lib.file_io import get_seq_content
from lib.six_frame_translation import six_frame_translate, extract_frame_aa, codon2aa


class SeqData():
    def __init__(self, filename):
        self._data = get_seq_content(filename)

    def get_data(self):
        return self._data

    def calculate_seq_length(self):
        return {header: len(seq) for header, seq in self._data.items()}

    def translate_content(self, frame):
        for header in self._data:
            self._data[header] = six_frame_translate(self._data[header], frame)


class BlastData():
    def __init__(self, filename):
        self._filename = filename
        self._data = defaultdict(lambda: '')
        self.read_blast_out()

    def read_blast_out(self):
        with open(self._filename, 'r') as f:
            for line in f:
                elements = line.split('\n')[0].split('\t')
                if self._data[elements[0]] == '':
                    self._data[elements[0]] = [elements[1:]]
                else:
                    self._data[elements[0]].append(elements[1:])

    def get_data(self):
        return self._data

    def get_best_hit(self):
        for allele, stats in self._data.items():
            # Sort by bitscore
            sorted(self._data[allele], key=lambda stat: float(stat[10]), reverse=False)
            # Sort by identity
            sorted(self._data[allele], key=lambda stat: float(stat[1]), reverse=False)
            # Sort by alignment length
            sorted(self._data[allele], key=lambda stat: float(stat[2]), reverse=False)

        return {allele: stats[0] for allele, stats in self._data.items()}
