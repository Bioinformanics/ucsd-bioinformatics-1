import sys

lines = sys.stdin.read().splitlines()
dna = lines[0].strip()

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

print(' '.join(str(index) for index in get_minimum_skews(dna)))

def get_minimum_skews_with_unknown_error(dna):
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
    return min_skew_indices
