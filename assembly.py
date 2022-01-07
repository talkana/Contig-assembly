#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from os import path
from lecture_functions import DeBruijnGraph


def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('reads_path', type=str, help="Path to input file containing reads in fasta format")
    parser.add_argument('contigs_path', type=str,
                        help="Path to output file containing assembled contigs in fasta format")
    args = parser.parse_args()
    return args


def check_input(reads_path, contigs_path):
    if not path.isfile(reads_path):
        raise Exception(f"File {reads_path} doesn't exist.")
    if not any(SeqIO.parse(reads_path, "fasta")):
        raise Exception(f"It looks like input file at path {reads_path} is empty or not in fasta format")
    if path.exists(contigs_path):
        raise Exception(f"Output path {contigs_path} already exists. Refusing to overwrite.")


def correct_errors():
    pass


def refine():
    """ Remove remaining "islands", “tips” and “bubbles” so that contigs are more obvious"""
    pass
