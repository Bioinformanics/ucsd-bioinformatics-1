import sys

lines = sys.stdin.read().splitlines() # read in the input from STDIN
dna = str.strip(lines[0])
k = int(str.strip(lines[1]))

def get_most_freq_k_mer(dna, k):
    dna = str.upper(dna)
    length = len(dna)
    dict = {}
    for index in range(0, length - k + 1):
        k_mer = dna[index: index + k]
        if dict.__contains__(k_mer):
            dict[k_mer] += 1
        else:
            dict[k_mer] = 1

    max_count = 0
    most_freq_k_mers = []
    for k_mer, count in iter(dict.items()):
        if count > max_count:
            most_freq_k_mers = [k_mer]
            max_count = count
        elif count == max_count:
            most_freq_k_mers.append(k_mer)

    return most_freq_k_mers

print(' '.join(get_most_freq_k_mer(dna, k)))