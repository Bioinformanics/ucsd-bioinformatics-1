# coding=utf-8
import sys
from collections import OrderedDict

"""
    PatternCount(Text, Pattern)
        count ← 0
        for i ← 0 to |Text| − |Pattern|
            if Text(i, |Pattern|) = Pattern
                count ← count + 1
        return count
"""
def count_pattern(dna, pattern):
    dna = dna.upper()
    pattern = pattern.upper()
    count = 0
    start = 0
    while True:
        start = dna.find(pattern, start) + 1
        if start > 0:
            count += 1
        else:
            return count


"""
    FrequentWords(Text, k)
        FrequentPatterns ← an empty set
        for i ← 0 to |Text| − k
            Pattern ← the k-mer Text(i, k)
            Count(i) ← PatternCount(Text, Pattern)
        maxCount ← maximum value in array Count
        for i ← 0 to |Text| − k
            if Count(i) = maxCount
                add Text(i, k) to FrequentPatterns
        return FrequentPatterns
"""
def get_most_freq_n_mer(dna, n):
    dna = str.upper(dna)
    length = str.__len__(dna)
    dict = {}
    for index in range(0, length - n):
        n_mer = dna[index: index + n]
        if dict.__contains__(n_mer):
            dict[n_mer] += 1
        else:
            dict[n_mer] = 1

    max_count = 0
    most_freq_n_mers = []
    for n_mer, count in dict.iteritems():
        if count > max_count:
            most_freq_n_mers = [n_mer]
            max_count = count
        elif count == max_count:
            most_freq_n_mers.append(n_mer)

    print 'Most Frequent n-mer Count: %d' % max_count
    for n_mer in most_freq_n_mers:
        print n_mer
    return most_freq_n_mers, max_count


def get_reverse_complement(dna):
    dna = str.upper(dna)
    complement = ''
    for base in dna:
        if base == 'A':
            complement += 'T'
        elif base == 'T':
            complement += 'A'
        elif base == 'C':
            complement += 'G'
        elif base == 'G':
            complement += 'C'
        else:
            raise Exception('Invalid DNA base.')
    return complement[::-1]


'''
CODE CHALLENGE: Solve the Pattern Matching Problem.
     Input: Two strings, Pattern and Genome.
     Output: A collection of space-separated integers specifying all starting positions where Pattern appears
     as a substring of Genome.
'''
def match_pattern(dna, pattern, match_rc=False):
    dna = dna.upper()
    pattern = pattern.upper()
    pattern_rc = get_reverse_complement(pattern)
    locations = []
    for index in range(len(dna) - len(pattern)):
        s = dna[index:index + len(pattern)]
        if s == pattern or (match_rc and s == pattern_rc):
            locations.append(str(index))
    print ", ".join(locations)
    return locations


'''
Clump Finding Problem: Find patterns forming clumps in a string.
    Input: A string Genome, and integers k, L, and t.
    Output: All distinct k-mers forming (L, t)-clumps in Genome.
Definition of Clump:
    Given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome
    if there is an interval of Genome of length L in which this k-mer appears at least t times.
    (This definition assumes that the k-mer completely fits within the interval.)
'''
def find_clump(dna, k, t, L):
    s = dna[0: L]
    counts = _count_frequencies(s, k)
    frequent_set = set()
    for pattern, count in counts.iteritems():
        if count >= t:
            frequent_set.add(pattern)

    for index in range(1, len(dna)-L+1):
        begin = dna[index - 1: index - 1 + k]
        counts[begin] -= 1
        end = dna[index + L - k: index + L]
        increase_count(counts, end)
        for pattern, count in counts.iteritems():
            if count >= t:
                frequent_set.add(pattern)

    frequency_array = sorted(list(frequent_set))
    return frequency_array


def _count_frequencies(dna, k):
    dict = {}
    for index in range(len(dna)-k+1):
        pattern = dna[index: index + k]
        increase_count(dict, pattern)
    return dict


def increase_count(dict, pattern):
    if dict.__contains__(pattern):
        dict[pattern] += 1
    else:
        dict[pattern] = 1
