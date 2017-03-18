from Week1.week1_utility import *

def find_clump_ecoli():
    with open('Datasets/E_coli.txt', 'r') as datafile:
        dna = datafile.readline()
    clumps = find_clump_efficient(dna, 9, 3, 500)
    print('Clumps of E-coli: %s' % str(clumps))
    print('Clump #: %s' % str(len(clumps)))


find_clump_ecoli()
#print str(pattern_to_number('ATGCAA'))
#print number_to_pattern(5437, 8)
#print ' '.join(map(lambda x: str(x), compute_frequencies('ACGCGGCTCTGAAA', 2)))

def frequency_array_quiz():
    with open('Datasets/dataset_2994_5.txt', 'r') as datafile:
        text = datafile.readlines()
        dna = str.strip(text[0])
        k = int(str.strip(text[1]))
    print(' '.join(map(lambda x: str(x), compute_frequencies(dna, k))))

def pattern_to_number_quiz():
    with open('Datasets/dataset_3010_2.txt', 'r') as datafile:
        pattern = str.strip(datafile.readline())
    print(pattern_to_number(pattern))

def number_to_pattern_quiz():
    with open('Datasets/dataset_3010_5.txt', 'r') as datafile:
        text = datafile.readlines()
        number = int(str.strip(text[0]))
        k = int(str.strip(text[1]))
    print(number_to_pattern(number, k))

#frequency_array_quiz()
#pattern_to_number_quiz()
#number_to_pattern_quiz()
