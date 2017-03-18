from Week1.week1_utility import get_reverse_complement
import itertools

def skew(dna):
    values = [0]
    for nucleotide in dna:
        value = values[-1]
        if nucleotide == 'C':
            value -= 1
        elif nucleotide == 'G':
            value += 1
        values.append(value)
    return values


def get_minimum_skews(dna):
    """
    Find a position in a genome where the skew diagram attains a minimum.
    :param dna: A DNA string Genome.
    :return: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    skew_values = skew(dna)
    min_value_indices = []
    min_value = 1000000
    for index in range(len(dna)+1):
        value = skew_values[index]
        if value < min_value:
            min_value_indices = [index]
            min_value = value
        elif value == min_value:
            min_value_indices.append(index)
    return min_value_indices
    """
    skew = 0
    min_skew = 10
    min_skew_indices = []
    for index in range(len(dna)):
        nucleotide = dna[index]
        if nucleotide == 'C':
            skew -= 1
            if skew < min_skew:
                min_skew = skew
                min_skew_indices = [index+1]
            elif skew == min_skew:
                min_skew_indices.append(index+1)
        elif nucleotide == 'G':
            skew += 1
    return min_skew_indices"""


def get_maximum_skews(dna):
    """
    Find a position in a genome where the skew diagram attains a maximum.
    :param dna: A DNA string Genome.
    :return: All integer(s) i maxmizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    skew_values = skew(dna)
    max_value_indices = []
    max_value = -1000000
    for index in range(len(dna)+1):
        value = skew_values[index]
        if value > max_value:
            max_value_indices = [index]
            max_value = value
        elif value == max_value:
            max_value_indices.append(index)
    return max_value_indices


def get_hamming_distance(dna1, dna2):
    """
    Hamming Distance Problem: Compute the Hamming distance between two strings.
        Input: Two strings of equal length.
    :return: The Hamming distance between these strings.
    """
    if len(dna1) != len(dna2):
        raise Exception("The two DNA strings must have equal length.")

    hamming_distance = 0
    for index in range(len(dna1)):
        if dna1[index] != dna2[index]:
            hamming_distance += 1
    return hamming_distance


def approximate_pattern_matching(pattern, dna, mismatches):
    """
    Find all approximate occurrences of a pattern in a string.
    :param pattern: the k-mer pattern string
    :param dna: the DNA string where pattern is searched against
    :param mismatches: the maximum hamming distance allowed for a match
    :return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    """
    pattern_len = len(pattern)
    matches = []
    for index in range(len(dna)-pattern_len+1):
        sub_dna = dna[index:index+pattern_len]
        hamming_distance = 0
        is_match = True
        for i in range(pattern_len):
            if pattern[i] != sub_dna[i]:
                hamming_distance += 1
                if hamming_distance > mismatches:
                    is_match = False
                    break
        if is_match:
            matches.append(index)
    return matches


def approximate_pattern_count(pattern, dna, mismatches):
    """
    Count occurrences of k-mers, possibly with mismatches.
    :param pattern: the k-mer pattern string
    :param dna: the DNA string where pattern is searched against
    :param mismatches: the maximum hamming distance allowed for a match
    :return: Count occurrence of Pattern appears as a substring of a DNA with at most d mismatches.
    """
    return len(approximate_pattern_matching(pattern, dna, mismatches))


def frequent_words_with_mismatches(dna, k, d):
    """
    Find the most frequent k-mers with mismatches in a string.
     Input: A string Text as well as integers k and d. (You may assume k <= 12 and d <= 3.)
    :param dna: the DNA string
    :param k: length of k-mer
    :param d: the maximum hamming distance allowed for a match
    :return: All most frequent k-mers with up to d mismatches in Text.
    """
    array_size = 4 ** k
    dna_length = len(dna)
    close = [0] * array_size
    max_count = 0
    for i in range(dna_length - k + 1):
        neighborhood = get_neighbours(dna[i:i + k], d)
        for pattern in neighborhood:
            index = _pattern_to_number(pattern)
            count = close[index] + 1
            close[index] = count
            if count > max_count:
                max_count = count

    max_indices = []
    for i in range(array_size):
        if close[i] == max_count:
            max_indices.append(i)

    return [_number_to_pattern(index, k) for index in max_indices]


def _pattern_to_number(pattern):
    number = 0
    for base in pattern:
        number *= 4
        number += _base_ids[base]
    return number


def _number_to_pattern(number, size):
    reverse_array = []
    for i in range(size):
        residue = number % 4
        number //= 4
        reverse_array.append(_bases[residue])
    return ''.join(reverse_array[::-1])


_bases = ['A', 'C', 'G', 'T']
_base_ids = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def get_neighbours(pattern, d):
    """
    For a given k-mer Pattern, gets all d-neighborhood with Hamming distance not exceeding d from this k-mer.
    :param pattern: the k-mer pattern
    :param d: maximum Hamming distance allowed
    :return: all d-neighborhood, including the given pattern
    """
    neighbours = set()
    _find_all_neighbours(pattern, '', d, neighbours)
    return neighbours

def _find_all_neighbours(pattern, prefix, d, neighbours):
    if not pattern: return
    sub_pattern = pattern[1:]
    for base in _bases:
        delta_d = 0 if base == pattern[0] else 1
        new_d = d - delta_d
        if new_d < 0:
            continue
        new_prefix = prefix + base
        neighbours.add(new_prefix + sub_pattern)
        if sub_pattern:
            _find_all_neighbours(sub_pattern, new_prefix, d - delta_d, neighbours)


def frequent_words_with_mismatches_and_reverse_complements(dna, k, d):
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.
     Input: A string Text as well as integers k and d. (You may assume k <= 12 and d <= 3.)
    :param dna: the DNA string
    :param k: length of k-mer
    :param d: the maximum hamming distance allowed for a match
    :return: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Pattern reverse complement)
             over all possible k-mers.
    """
    array_size = 4 ** k
    dna_length = len(dna)
    close = [0] * array_size
    max_count = 0
    max_indices = []
    for i in range(dna_length - k + 1):
        pattern = dna[i:i + k]
        neighbors = get_neighbours(pattern, d)
        rc = get_reverse_complement(pattern)
        rc_neighbours = get_neighbours(rc, d)
        for p in itertools.chain(neighbors, rc_neighbours):
            index = _pattern_to_number(p)
            close[index] += 1
            count = close[index]
            if count > max_count:
                max_count = count
                max_indices = [index]
            elif count == max_count:
                max_indices.append(index)

    return [_number_to_pattern(index, k) for index in max_indices]


