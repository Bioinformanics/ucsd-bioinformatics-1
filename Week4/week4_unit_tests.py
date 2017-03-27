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
            dnas = []
            for i in range(t):
                dnas.append(datafile.readline().strip())
            datafile.readline() # "Output"
            expected_motifs = []
            for i in range(t):
                expected_motifs.append(datafile.readline().strip())

        motifs = randomized_motif_search(dnas, k, t)
        self.assertTrue(AreStringListsEqual(expected_motifs, motifs))

    def test_extra_dataset(self):
        self._test('Datasets/RandomizedMotifSearch/extra.txt')

    def test_sample_dataset(self):
        self._test('Datasets/RandomizedMotifSearch/sample.txt')

    def test_dataset_1(self):
        self._test('Datasets/RandomizedMotifSearch/data01.txt')

    def test_dataset_2(self):
        self._test('Datasets/RandomizedMotifSearch/data02.txt')

