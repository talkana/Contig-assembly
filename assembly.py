#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from os import path
from lecture_functions import kmerHist, correct1mm, DeBruijnGraph
from math import ceil
from memory_test import display_top
import tracemalloc
from collections import Counter


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('reads_path', type=str,
                        help="Path to input file containing reads in fasta format")
    parser.add_argument('contigs_path', type=str,
                        help="Path to output file containing assembled contigs in fasta format")
    args = parser.parse_args()
    return args.reads_path, args.contigs_path


def check_input(reads_path, contigs_path):
    if not path.isfile(reads_path):
        raise Exception(f"File {reads_path} doesn't exist.")
    if not any(SeqIO.parse(reads_path, "fasta")):
        raise Exception(f"It looks like input file at path {reads_path} is empty or not in fasta format")
    if path.exists(contigs_path):
        raise Exception(f"Output path {contigs_path} already exists. Refusing to overwrite.")


def get_reads(reads_path):
    record = SeqIO.parse(reads_path, "fasta")
    reads = [seq.seq for seq in record]  # todo: change to iterator to improve memory
    return reads


def select_parameters():
    """Returns kmer size and frequency threshold for kmer correction """
    avg_read_length = 80
    avg_read_number = 1000
    k = 7
    expected_kmer_occ = (avg_read_length - k + 1) * avg_read_number / 4 ** k
    freq_thr = ceil(expected_kmer_occ / 2)  # todo: check if optimal (plot histogram?)
    print(
        f"Expected number of kmer occurences is {expected_kmer_occ}. Will try to correct kmers with <= {freq_thr} occurences")
    return k, freq_thr


def correct_reads(reads, k, freq_threshold):
    kmerhist = kmerHist(reads, k)
    alphabet = ["A", "T", "G", "C"]
    corrected_reads = []  # todo: change to iterator to improve memory
    for read in reads:
        for letter in read:
            assert letter in alphabet
        read = correct1mm(read, k, kmerhist, alphabet, freq_threshold)
        corrected_reads.append(read)
    return corrected_reads


def refine(graph):
    """ optional """
    """ Remove remaining "islands", “tips” and “bubbles” so that contigs are more obvious"""
    pass


def extract_contigs(graph):
    """ Extracts contigs from de Bruijn graphs using the greedy approach (as in the greedy SCS algorithm)."""
    pass


def main():
    tracemalloc.start() # memory check
    counts = Counter() # memory check
    reads_path, contigs_path = parse_options()
    check_input(reads_path, contigs_path)
    reads = get_reads(reads_path)
    k, freq_thr = select_parameters()
    corrected_reads = correct_reads(reads, k, freq_thr)
    DB_graph = DeBruijnGraph(corrected_reads, k)
    snapshot = tracemalloc.take_snapshot() # memory check
    display_top(snapshot) # memory check


if __name__ == "__main__":
    main()
