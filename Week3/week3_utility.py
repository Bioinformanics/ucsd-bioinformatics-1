# coding=utf-8
import itertools
from Week2.week2_utility import get_neighbours, get_hamming_distance
from math import log
from utilities import reduce
from utilities import ConvertBaseToIndex, ConvertIndexToBase


def _find_motifs(dna, k, d):
    """
    Find all (k, d)-motifs in a dna string
    :param dna: the dna string
    :param k: k-mer
    :param d: the maximally allowed distance for a match
    :return: All (k, d)-motifs in Dna
    """
    motifs = set()
    for index in range(len(dna) - k + 1):
        motif = dna[index: index + k]
        neighbours = get_neighbours(motif, d)
        motifs.add(motif)
        for neighbour in neighbours:
            motifs.add(neighbour)
    return motifs


def motif_enumeration(dnas, k, d):
    """
    Find all (k, d)-motifs in a list of dna strings. This is a brute force algorithm to find motifs.
    :param dnas: a list of DNA strings
    :param k: k-mer
    :param d: the maximally allowed distance for a match
    :return: All (k, d)-motifs in Dna
    """
    motifs = []
    for dna in dnas:
        motifs.append(_find_motifs(dna, k, d))
    from operator import and_
    return reduce(and_, motifs)


def get_minimum_hamming_distance(pattern, dna):
    """
    Hamming Distance Problem: Compute the minimum Hamming distance between a pattern and a DNA string.
        Input: a pattern string and a dna string (which is usually longer than pattern)
    :param pattern: a DNA pattern string
    :param dna: A DNA string, usually longer than the pattern
    :return: The minimum Hamming distance.
    """
    k = len(pattern)
    min_distance = k
    for index in range(len(dna) - k + 1):
        sub_dna = dna[index:index + k]
        distance = get_hamming_distance(pattern, sub_dna)
        if distance < min_distance:
            min_distance = distance
    return min_distance


def find_median_string(dnas, k):
    """
    Median String Problem: Find a median string.
    :param dnas: A collection of strings Dna
    :param k: length of k-mer
    :return: All k-mer Patterns that minimizes d(Pattern, Dna).
    """
    ms = []
    min_d = len(dnas) * k
    for l in itertools.product('ACGT', repeat=k):
        k_mer = ''.join(l)
        d = reduce(lambda a, b: a + b, map(lambda dna: get_minimum_hamming_distance(k_mer, dna), dnas))
        if d < min_d:
            min_d = d
            ms = [k_mer]
        elif d == min_d:
            ms.append(k_mer)
    return ms


def get_profile_most_probable_k_mer(dna, profile_matrix):
    """
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
    (https://stepik.org/lesson/Greedy-Motif-Search-159/step/3?course=Stepic-Interactive-Text-for-Week-3&unit=8217)
    :param dna: a DNA string
    :param profile_matrix: a 4 × k matrix Profile
    :return: A Profile-most probable k-mer in Text.
    """
    k = len(profile_matrix[0])
    max_probability = 0
    max_k_mer = None
    for index in range(len(dna) - k + 1):
        k_mer = dna[index: index + k]
        probability = 1
        for pos, base in enumerate(k_mer):
            probability *= profile_matrix[ConvertBaseToIndex(base)][pos]
        if probability > max_probability:
            max_probability = probability
            max_k_mer = k_mer
    return max_k_mer

# Greedy motif search: https://github.com/jschendel/Rosalind/blob/master/Textbook_03D.py

def profile_most_probable_kmer(dna, k, profile):
    '''Returns the profile most probable k-mer for the given input data.'''
    # A dictionary relating nucleotides to their position within the profile.
    nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}

    # Initialize the maximum probabily.
    max_probability = -1

    # Compute the probability of the each k-mer, store it if it's currently a maximum.
    for i in range(len(dna)-k+1):
        # Get the current probability.
        current_probability = 1
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_probability *= profile[j][nuc_loc[nucleotide]]

        # Check for a maximum.
        if current_probability > max_probability:
            max_probability = current_probability
            most_probable = dna[i:i+k]

    return most_probable

def score(motifs):
    '''Returns the score of the given list of motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    max_count = sum([max([c.count(nucleotide) for nucleotide in 'ACGT']) for c in columns])
    return len(motifs[0])*len(motifs) - max_count


def profile(motifs):
    '''Returns the profile of the dna list motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    return [[float(col.count(nuc)) / float(len(col)) for nuc in 'ACGT'] for col in columns]


def greedy_motif_search(dna_list, k):
    '''Runs the Greedy Motif Search algorithm and returns the best motif.'''
    # Initialize the best score as a score higher than the highest possible score.
    t = len(dna_list)
    best_score = t*k
    best_motifs = None

    # Run the greedy motif search.
    for i in range(len(dna_list[0])-k+1):
        # Initialize the motifs as each k-mer from the first dna sequence.
        motifs = [dna_list[0][i:i+k]]

        # Find the most probable k-mer in the next string.
        for j in range(1, t):
            current_profile = profile(motifs)
            motifs.append(profile_most_probable_kmer(dna_list[j], k, current_profile))

        # Check to see if we have a new best scoring list of motifs.
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs

    return best_motifs


def profile_laplace(motifs):
    '''Returns the profile of the dna list motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    return [[float(col.count(nuc) + 1) / float(len(col) + 4) for nuc in 'ACGT'] for col in columns]


def greedy_motif_search_laplace(dna_list, k):
    """
    Runs the Greedy Motif Search algorithm and returns the best motif,
    with applications of Laplace's Rule of Succession
    """
    # Initialize the best score as a score higher than the highest possible score.
    t = len(dna_list)
    best_score = t*k
    best_motifs = None

    # Run the greedy motif search.
    for i in range(len(dna_list[0])-k+1):
        # Initialize the motifs as each k-mer from the first dna sequence.
        motifs = [dna_list[0][i:i+k]]
        current_profile = profile_laplace(motifs)

        # Find the most probable k-mer in the next string.
        for j in range(1, t):
            motifs.append(profile_most_probable_kmer(dna_list[j], k, current_profile))
            current_profile = profile_laplace(motifs)

        # Check to see if we have a new best scoring list of motifs.
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs

    return best_motifs

