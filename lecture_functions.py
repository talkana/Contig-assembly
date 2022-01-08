import matplotlib.pyplot as plt
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
        self.weights = {} # maps graph edge (from, to) to edge weight
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st, k):
                self.k1mers.update([km1L, km1R])
                if not (km1L, km1R) in self.weights.keys():
                    self.weights[(km1L, km1R)] = 1
                else:
                    self.weights[(km1L, km1R)] += 1


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
    frequencies = list(kmerhist.values())
    plt.hist(frequencies)
    plt.show()
    return kmerhist


def correct1mm(read, k, kmerhist, alpha, thresh):
    """ Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. """
    # Iterate over k-mers in read
    for i in range(len(read) - (k - 1)):
        kmer = read[i:i + k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i + k:]
                    break
    # Return possibly-corrected read
    return read


