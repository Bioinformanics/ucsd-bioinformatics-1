import sys

lines = sys.stdin.read().splitlines() # read in the input from STDIN
dna = lines[0].strip()
k = int(lines[1].strip())

_base_to_number_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def pattern_to_number(pattern):
    number = 0
    for base in pattern:
        number = number * 4 + _base_to_number_dict[base]
    return number

def compute_frequencies(dna, k):
    frequency_array = [0] * 4**k
    for index in range(0, len(dna) - k + 1):
        pattern = dna[index:index + k]
        array_index = pattern_to_number(pattern)
        frequency_array[array_index] += 1
    return frequency_array


print(" ".join(map(lambda x: str(x), compute_frequencies(dna, k))))
