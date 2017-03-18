import sys

lines = sys.stdin.read().splitlines() # read in the input from STDIN
index = int(lines[0].strip())
k = int(lines[1].strip())

_number_to_base_table = ['A', 'C', 'G', 'T']

def number_to_pattern(number, k):
    pattern = []
    while k>0:
        remainder = number % 4
        number //= 4
        pattern.insert(0, _number_to_base_table[remainder])
        k -= 1
    return ''.join(pattern)

print(number_to_pattern(index, k))