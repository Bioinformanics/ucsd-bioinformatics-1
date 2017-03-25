# coding=utf-8

from Week2.week2_utility import get_hamming_distance
from utilities import ConvertBaseToIndex, ConvertIndexToBase

class Profile:
    """
    Example: https://stepik.org/lesson/Greedy-Motif-Search-159/step/2?course=Stepic-Interactive-Text-for-Week-3&unit=8217
    """
    def __init__(self, motifs):
        self._construct_pr(motifs)
        self._construct_consensus()
        self._calculate_score(motifs)

    def _construct_pr(self, motifs):
        motif_length = len(motifs[0])
        self._pr = []
        for j in range(4):
            self._pr.append([0] * motif_length)
        for i in range(motif_length):
            counts = [0] * 4
            for motif in motifs:
                base_index = ConvertBaseToIndex(motif[i])
                counts[base_index] += 1
            total = sum(counts)
            for base_id in range(4):
                self._pr[base_id][i] = counts[base_id] / total

    def _construct_consensus(self):
        motif_length = len(self._pr[0])
        self._consensus = []
        for i in range(motif_length):
            max_base_index = -1
            max_pr = 0
            for base_index in range(4):
                if self._pr[base_index][i] > max_pr:
                    max_pr = self._pr[base_index][i]
                    max_base_index = base_index
            self._consensus.append(ConvertIndexToBase(max_base_index))

        #selfs.entropy = sum(map(lambda p: p * log(p, 2) if p > 0 else 0, self._probabilities))
        #if self.entropy < 0:
        #    self.entropy *= -1

    def _calculate_score(self, motifs):
        self.score = 0
        for motif in motifs:
            self.score += get_hamming_distance(motif, self._consensus)
        return self.score


def score(motifs):
    return Profile(motifs).score


"""
    GreedyMotifSearch(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
"""
def greedy_motif_search(dna_list, k):
    """
    Find a list of best motif (k-mer) which minimizes the overall scoring.
    If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
    :param dna_list: a list of DNA
    :param k: length of k-mer
    :return: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t).
    """
    t = len(dna_list)
    dna_length = len(dna_list[0])
    best_motifs = [dna[0:k] for dna in dna_list]
    best_motifs_score = score(best_motifs)
    for index in range(0, dna_length - k + 1):
        motifs = [dna[index:index+k] for dna in dna_list]
        motifs_score = score(motifs)
        if motifs_score < best_motifs_score:
            best_motifs_score = motifs_score
            best_motifs = motifs
    return best_motifs