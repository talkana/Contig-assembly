#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from os import path


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
    reads = [str(seq.seq) for seq in record]
    return reads


def select_parameters():
    """Return recommended kmer size and frequency threshold for kmer correction """
    k = 18
    freq_thr = 1
    return k, freq_thr


def neighbors1mm(kmer, alpha):
    """ Generate all neighbors at Hamming distance 1 from kmer """
    neighbors = []
    for j in range(len(kmer) - 1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j + 1:])
    return neighbors


def kmerHist(reads, k):
    """ Return k-mer histogram and average k-mer occurrences """
    kmerhist = {}
    for read in reads:
        for kmer in [read[i:i + k] for i in range(len(read) - (k - 1))]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist


def correct1mm(read, k, kmerhist, alpha, thresh):
    """ Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. """
    for i in range(len(read) - (k - 1)):
        kmer = read[i:i + k]
        if kmerhist.get(kmer, 0) <= thresh:
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    read = read[:i] + newkmer + read[i + k:]
                    break
    return read


def correct_reads(reads, k, freq_threshold):
    kmerhist = kmerHist(reads, k)
    alphabet = ["A", "T", "G", "C"]
    corrected_reads = []
    for read in reads:
        read = correct1mm(read, k, kmerhist, alphabet, freq_threshold)
        corrected_reads.append(read)
    return corrected_reads


class DeBruijnGraph:
    """ De Bruijn directed weighted graph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(st, k):
        """ Chop string into k-mers of given length """
        for i in range(len(st) - (k - 1)):
            yield st[i:i + k], st[i:i + k - 1], st[i + 1:i + k]

    def __init__(self, strIter, k):
        """ Build de Bruijn weighted graph given string iterator and k-mer
            length k """
        self.k1mers = set()
        self.weights = {}  # maps graph edge (from, to) to edge weight
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st, k):
                self.k1mers.update([km1L, km1R])
                if not (km1L, km1R) in self.weights.keys():
                    self.weights[(km1L, km1R)] = 1
                else:
                    self.weights[(km1L, km1R)] += 1


def get_contigs_greedy(DB_graph, k):
    """Get contigs from DeBruijn graph using greedy approach: iteratively merge pairs of nodes,
     starting with the pair that is connected by the edge with highest weight"""

    contigs = list(DB_graph.k1mers)
    weights = DB_graph.weights  # map from (Node1, Node2) to edge weight
    weights = sorted([(item[0][0], item[0][1], item[1]) for item in weights.items()],
                     key=lambda x: x[2])  # sorted by weight

    while weights:
        best_edge = weights.pop()
        left_seq = best_edge[0]
        right_seq = best_edge[1]
        new_contig = left_seq + right_seq[k - 2:]
        contigs.remove(left_seq)
        contigs.remove(right_seq)
        contigs.append(new_contig)
        new_weights = []
        for edge in weights:  # update edges
            el = edge[0]
            er = edge[1]
            if el == er:
                pass
            elif el == right_seq and er == left_seq:  # would create a cycle, delete
                pass
            elif el == right_seq:  # edges starting in the 2nd node should start in the merged node
                new_edge = (new_contig, edge[1], edge[2])
                new_weights.append(new_edge)
            elif er == left_seq:  # edges pointing to the 1st node should point to the merged node
                new_edge = (edge[0], new_contig, edge[2])
                new_weights.append(new_edge)
            elif el != left_seq and er != right_seq:  # delete edges starting at 1st node or pointing to 2nd node
                new_weights.append(edge)
        weights = new_weights
    return contigs


def filter_contigs(contigs, minlen):
    filtered = []
    for contig in contigs:
        if len(contig) >= minlen:
            filtered.append(contig)
    return filtered


def save_contigs(contigs, filename):
    contigfile = open(filename, "w")
    for i in range(len(contigs)):
        if len(contigs[i]) > 300:
            contigfile.write(f">Contig{i}\n")
            contigfile.write(f"{contigs[i]}\n")
    contigfile.close()


def main():
    reads_path, contigs_path = parse_options()
    check_input(reads_path, contigs_path)
    reads = get_reads(reads_path)
    k, freq_thr = select_parameters()
    corrected_reads = correct_reads(reads, k, freq_thr)
    DB_graph = DeBruijnGraph(corrected_reads, k)
    contigs = get_contigs_greedy(DB_graph, k)
    long_contigs = filter_contigs(contigs, 300)
    save_contigs(long_contigs, contigs_path)


if __name__ == "__main__":
    main()
