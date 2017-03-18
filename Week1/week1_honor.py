from week1_utility import *

_base_to_number_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_number_to_base_table = ['A', 'C', 'G', 'T']

def pattern_to_number(pattern):
    number = 0
    for base in pattern:
        number = number * 4 + _base_to_number_dict[base]
    return number

def number_to_pattern(number, k):
    pattern = []
    while k>0:
        remainder = number % 4
        number /= 4
        pattern.insert(0, _number_to_base_table[remainder])
        k -= 1
    return ''.join(pattern)

def compute_frequencies(text, k):
    frequency_array_size = 4**k
    frequency_array = [0] * frequency_array_size
    for index in range(0, len(text) - k + 1):
        pattern = text[index:index+k]
        array_index = pattern_to_number(pattern)
        frequency_array[array_index] += 1
    return frequency_array

def find_clump_ecoli():
    with open('E_coli.txt', 'r') as datafile:
        dna = datafile.readline()
    clumps = find_clump(dna, 9, 3, 500)
    print 'Clumps of E-coli: %s' % str(clumps)
    print 'Clump #: %s' % str(len(clumps))


#find_clump_ecoli()
#print str(pattern_to_number('ATGCAA'))
#print number_to_pattern(5437, 8)
#print ' '.join(map(lambda x: str(x), compute_frequencies('ACGCGGCTCTGAAA', 2)))

def frequency_array_quiz():
    with open('Datasets/dataset_2994_5.txt', 'r') as datafile:
        text = datafile.readlines()
        dna = str.strip(text[0])
        k = int(str.strip(text[1]))
    print ' '.join(map(lambda x: str(x), compute_frequencies(dna, k)))

def pattern_to_number_quiz():
    with open('Datasets/dataset_3010_2.txt', 'r') as datafile:
        pattern = str.strip(datafile.readline())
    print pattern_to_number(pattern)

def number_to_pattern_quiz():
    with open('Datasets/dataset_3010_5.txt', 'r') as datafile:
        text = datafile.readlines()
        number = int(str.strip(text[0]))
        k = int(str.strip(text[1]))
    print number_to_pattern(number, k)

#frequency_array_quiz()
#pattern_to_number_quiz()
number_to_pattern_quiz()
