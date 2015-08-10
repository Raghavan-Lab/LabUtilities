#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Splits a multi-sequence fasta file into multiple files containing n number of reads.

   Examples:
   Split a very large fasta file into smaller files containing 10,000 reads
     python split_fasta.py -f large_amino_acid_file.faa -n 10000


   NOTE: This script places the output files in the current directory.  It's recommended that
   you create an empty directory to work in.

"""


import argparse
from collections import namedtuple

from Bio import SeqIO


def get_options():
    """
    Returns command line arguments as an Argparse object
    :return: command line options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help='Number of reads per output file', default=1000)
    parser.add_argument('-f', help='Fasta file to split')
    return parser.parse_args()


def save_fasta_slice(slice, start_record, end_record):
    """
    Saves a list of namedtuples as a fasta file

    :param slice: list of reads to save
    :type : list

    :param start_record: first record in slice
     :type : int

    :param end_record: last record in slice
    :type : int

    """
    output_file = 'sequences_{}_to_{}.fasta'.format(start_record, end_record)
    with open(output_file, 'w') as fasta_slice:
        for sequence in slice:
            # output is in fasta format
            fasta_slice.write('{id}\n{seq}\n\n'.format(id=sequence.seq_id,
                                                       seq=sequence.seq))


def main():
    """
    Reads a fasta file into a list of tuples (sequence id, nucleotide/amino acid sequence) and
    saves n sized batches into output files.
    """
    # get command line arguments
    options = get_options()

    # Keep track of the read numbers for each chunk
    first_seq = 1
    last_seq = 0

    sequence_buffer = []
    FastaRead = namedtuple('FastaRead', ['seq_id', 'seq'])
    for record in SeqIO.parse(options.f, "fasta"):
        last_seq += 1
        sequence_buffer.append(FastaRead('>{}'.format(record.description), record.seq))
        # Save a file if we have enough reads
        if len(sequence_buffer) >= int(options.n):
            # flush the buffer
            save_fasta_slice(sequence_buffer, first_seq, last_seq)
            print('\tsaved records: {} to {}'.format(first_seq, last_seq))
            # Reset buffer and starting record number
            sequence_buffer = []
            first_seq = last_seq + 1
    # save remaining bits
    save_fasta_slice(sequence_buffer, first_seq, last_seq)


if __name__ == '__main__':
    print('Splitting Fasta file...')
    main()
    print('Job Complete')
