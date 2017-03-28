import unittest
from Week4.week4_utility import *
from utilities import AreStringListsEqual

class TestRandomizedMotifSearch(unittest.TestCase):
    def _test(self, datafile_name):
        with open(datafile_name, 'r') as datafile:
            datafile.readline() # "Input"
            args = datafile.readline().strip().split(" ")
            k = int(args[0])
            t = int(args[1])
            dna_list = []
            for i in range(t):
                dna_list.append(datafile.readline().strip())
            datafile.readline() # "Output"
            expected_motifs = []
            for i in range(t):
                expected_motifs.append(datafile.readline().strip())

        motifs = randomized_motif_search(dna_list, k, t)
        passed = AreStringListsEqual(expected_motifs, motifs)
        if not passed:
            print("Expected: " + ' '.join(expected_motifs))
            print("Actual: " + ' '.join(motifs))
        self.assertTrue(passed)

    def test_extra_dataset(self):
        self._test('Datasets/RandomizedMotifSearch/extra.txt')

    def test_sample_dataset(self):
        self._test('Datasets/RandomizedMotifSearch/sample.txt')

    def test_dataset_1(self):
        self._test('Datasets/RandomizedMotifSearch/data01.txt')

    def test_dataset_2(self):
        self._test('Datasets/RandomizedMotifSearch/data02.txt')

