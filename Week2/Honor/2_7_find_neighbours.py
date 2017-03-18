import sys

lines = sys.stdin.read().splitlines()
pattern = lines[0].strip()
d = int(lines[1].strip())

_bases = ['A', 'C', 'G', 'T']

def get_neighbours(pattern, d):
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

print('\n'.join(get_neighbours(pattern, d)))
