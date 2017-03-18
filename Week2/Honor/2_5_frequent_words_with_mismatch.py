import sys

lines = sys.stdin.read().splitlines()
dna = lines[0].strip()
args = lines[1].strip().split(' ')
k = int(args[0])
d = int(args[1])

_bases = ['A', 'C', 'G', 'T']
_base_ids = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

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

def get_neighbours(pattern, d):
    neighbours = set()
    _find_all_neighbours(pattern, '', d, neighbours)
    return neighbours

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

def frequent_words_with_mismatches(dna, k, d):
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

print(' '.join(frequent_words_with_mismatches(dna, k, d)))
