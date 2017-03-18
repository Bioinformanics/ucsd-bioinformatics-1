import sys

lines = sys.stdin.read().splitlines() # read in the input from STDIN
dna = lines[0].strip()

_base_to_number_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def pattern_to_number(pattern):
    number = 0
    for base in pattern:
        number = number * 4 + _base_to_number_dict[base]
    return number

print(str(pattern_to_number(dna)))
