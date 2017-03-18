import sys

lines = sys.stdin.read().splitlines() # read in the input from STDIN
pattern = str.strip(lines[0])
dna = str.strip(lines[1])

def match_pattern(dna, pattern):
    dna = dna.upper()
    pattern = pattern.upper()
    locations = []
    for index in range(len(dna) - len(pattern) + 1):
        s = dna[index:index + len(pattern)]
        if s == pattern:
            locations.append(str(index))
    return locations

print(" ".join(match_pattern(dna, pattern)))