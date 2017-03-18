from week1_utility import *

def find_clump_ecoli():
    with open('E_coli.txt', 'r') as datafile:
        dna = datafile.readline()
    find_clump(dna, 9, 3, 500)

find_clump_ecoli()