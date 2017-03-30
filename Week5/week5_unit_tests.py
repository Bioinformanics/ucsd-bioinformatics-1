import unittest
from Week5.week5 import week5

class TestWeek5(unittest.TestCase):
    def test_week5(self):
        expected_motif_counts = []
        with open('expected_counts.txt', 'r') as file:
            lines = file.readlines()
            for line in lines:
                expected_motif_counts.append(line.strip().split(' '))

        motifs = week5()
        columns = [''.join(dna) for dna in zip(*motifs)]
        actual_counts = [[column.count(base) for base in 'ACGT'] for column in columns]

        for i in range(len(expected_motif_counts)):
            for j in range(4):
                equal = expected_motif_counts[i][j] == actual_counts[i],[j]
                if not equal:
                    k = 0
                self.assertTrue(equal)

