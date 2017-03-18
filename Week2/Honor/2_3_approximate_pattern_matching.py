import sys

lines = sys.stdin.read().splitlines()
pattern = lines[0].strip()
dna = lines[1].strip()
d = int(lines[2].strip())

def approximate_pattern_matching(pattern, dna, mismatches):
    pattern_len = len(pattern)
    matches = []
    for index in range(len(dna)-pattern_len+1):
        sub_dna = dna[index:index+pattern_len]
        hamming_distance = 0
        is_match = True
        for i in range(pattern_len):
            if pattern[i] != sub_dna[i]:
                hamming_distance += 1
                if hamming_distance > mismatches:
                    is_match = False
                    break
        if is_match:
            matches.append(index)
    return matches

print(' '.join(str(index) for index in approximate_pattern_matching(pattern, dna, d)))
