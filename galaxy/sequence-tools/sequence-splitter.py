#!/usr/bin/env python2
# coding: utf-8

import locale
import sys
import os
import argparse
from itertools import chain, islice
# BioPython
from Bio import SeqIO


def main():
    args = parse_arguments()
    # extract the basename and the file extension
    basename, file_extension = os.path.splitext(args.sequences)
    # extract the filename (with no extension)
    _, filename = os.path.split(basename)

    # split the sequences in chunks
    if (args.chunk_size):
        sequences_record = gen_sequence_record(args.sequences, args.format)
        chunks = gen_get_chunks_by_size(sequences_record, args.chunk_size)
    else:
        chunks = gen_get_chunks(args.sequences, args.format, args.nb_chunk)

    # Write the chunks in numbered files.
    write_chunks(chunks, args.output, filename, file_extension, args.format)


def gen_get_chunks(sequences_path, sequences_format, nb_chunk):

    # First record to count the sequences
    sequences_record_to_count = gen_sequence_record(
        sequences_path, sequences_format)
    # Get the number of sequences
    nb_sequences = get_nb_sequences(sequences_record_to_count)

    # Second record to that will be splitted
    sequences_to_split = gen_sequence_record(sequences_path, sequences_format)

    # Get the size of the chunks
    chunk_size = int(nb_sequences / nb_chunk) if nb_sequences > nb_chunk else 1 
    return gen_get_chunks_by_size(sequences_to_split, chunk_size)


def gen_get_chunks_by_size(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))


def gen_sequence_record(sequences_path, sequence_format):
    return SeqIO.parse(sequences_path, sequence_format)


def get_nb_sequences(sequences):
    count = 0
    for seq in sequences:
        count += 1
    return count


def write_chunks(iterable, dirname, filename, file_extension, sequence_format):
    for idx, chunk in enumerate(iterable):
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        with open(os.path.join(dirname, filename + '-chunk-' + str(idx+1) + file_extension), mode='w') as output_handle:
            SeqIO.write(chunk, output_handle, sequence_format)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Split fastq file")
    parser.add_argument("-s", "--sequences", type=str,
                        help="File that contains the sequences")

    parser.add_argument("-f", "--format", type=str,
                        help="File format (fastq, fasta)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-cs", "--chunk-size", type=int,
                       help="Size of the chunks. Cannot be used  ")
    group.add_argument("-nc", "--nb-chunk", help="Number of chunk", type=int)

    parser.add_argument("-o", "--output", type=str, default="./")

    return parser.parse_args()


if __name__ == '__main__':
    main()
