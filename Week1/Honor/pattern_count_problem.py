import sys

def pattern_count(text, pattern):
    count = 0
    for index in range(0, len(text)-len(pattern)+1):
        if text[index:index+len(pattern)] == pattern:
            count += 1
    return count

lines = sys.stdin.read().splitlines() # read in the input from STDIN
text = str.strip(lines[0])
pattern = str.strip(lines[1])

print(str(pattern_count(text, pattern)))