import sys

lines = sys.stdin.read().splitlines()
dna1 = lines[0].strip()
dna2 = lines[1].strip()

def get_hamming_distance(dna1, dna2):
    if len(dna1) != len(dna2):
        raise Exception("The two DNA strings must have equal length.")

    hamming_distance = 0
    for index in range(len(dna1)):
        if dna1[index] != dna2[index]:
            hamming_distance += 1
    return hamming_distance

print(str(get_hamming_distance(dna1, dna2)))