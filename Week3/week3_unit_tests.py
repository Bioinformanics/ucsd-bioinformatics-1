import unittest

from utilities import ConvertTextToMatrix, AreStringListsEqual
from Week3.week3_utility import *
from Week2.week2_unit_tests import ResultEqual

class TestMotifEnumeration(unittest.TestCase):
    def test_sample(self):
        k = 3
        d = 1
        dnas = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "ATA ATT GTT TTT"

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_1(self):
        k = 3
        d = 0
        dnas = ["ACGT", "ACGT", "ACGT"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "ACG CGT"

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_2(self):
        k = 3
        d = 1
        dnas = ["AAAAA", "AAAAA", "AAAAA"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "AAA AAC AAG AAT ACA AGA ATA CAA GAA TAA"

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_2_2(self):
        k = 3
        d = 0
        dnas = ["AAAAA", "AAAAA", "AAAAA"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "AAA"

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_3(self):
        k = 3
        d = 3
        dnas = ["AAAAA", "AAAAA", "AAAAA"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT"

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_4(self):
        k = 3
        d = 0
        dnas = ["AAAAA", "AAAAA", "AACAA"]
        motifs = motif_enumeration(dnas, k, d)
        expected = ""

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_5(self):
        k = 3
        d = 0
        dnas = ["AACAA", "AAAAA", "AAAAA"]
        motifs = motif_enumeration(dnas, k, d)
        expected = ""

        self.assertTrue(ResultEqual(expected, motifs))

    def test_dataset_6(self):
        k = 15
        d = 0
        dnas = ["atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg",
            "acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaaggggggga",
            "tgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccga",
            "gctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggaga",
            "tcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttatag",
            "gtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcacccaaattcagtgtgggcgagcgcaa",
            "cggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcat",
            "aacttgagttaaaaaaaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgta",
            "ttggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcataaaaaaaagggggggaccgaaagggaag",
            "ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga"]
        motifs = motif_enumeration(dnas, k, d)
        expected = "aaaaaaaaggggggg"

        self.assertTrue(ResultEqual(expected, motifs))

class TestFindMediumString(unittest.TestCase):
    def test_sample(self):
        dnas = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"]
        k = 3
        ms = find_median_string(dnas, k)
        expected = ["GAC"]
        self.assertEqual(expected, ms)

    def test_dataset_1(self):
        dnas = ["ACGT", "ACGT", "ACGT"]
        k = 3
        ms = find_median_string(dnas, k)
        expected = ["ACG", "CGT"]
        self.assertEqual(expected, ms)

    def test_dataset_2(self):
        dnas = "ATA ACA AGA AAT AAC".split(" ")
        k = 3
        ms = find_median_string(dnas, k)
        expected = ["AAA"]
        self.assertEqual(expected, ms)

    def test_dataset_3(self):
        dnas = "AAG AAT".split(" ")
        k = 3
        ms = find_median_string(dnas, k)
        expected = ["AAG", "AAT"]
        self.assertEqual(expected, ms)

    def test_extra_dataset(self):
        dnas = ["TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT",
                "CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA",
                "TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT",
                "TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA",
                "ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG",
                "TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA",
                "TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC",
                "GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA",
                "CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG",
                "CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG"]
        k = 6
        ms = find_median_string(dnas, k)
        expected = ["CGGCGA"]
        self.assertEqual(expected, ms)


class TestGetProfileMostProbableKmer(unittest.TestCase):
    def _test(self, datafile_name, expected):
        with open(datafile_name, 'r') as datafile:
            dna = datafile.readline().strip()
            k = int(datafile.readline().strip())
            matrix_text = []
            for loop in range(4):
                matrix_text.append(datafile.readline().strip())
        profile = ConvertTextToMatrix(matrix_text, float)
        k_mer = get_profile_most_probable_k_mer(dna, k, profile)
        self.assertEqual(expected, k_mer)

    def test_sample(self):
        self._test('Datasets/ProfileMostProbableKMer/sample.txt', "CCGAG")

    def test_dataset_1(self):
        self._test('Datasets/ProfileMostProbableKMer/data01.txt', "AGCAGCTT")

    def test_dataset_2(self):
        self._test('Datasets/ProfileMostProbableKMer/data02.txt', "AAGCAGAGTTTA")

    def test_debug_dataset(self):
        self._test('Datasets/ProfileMostProbableKMer/debug.txt', "TGTCGC")


class GreedyMotifSearch(unittest.TestCase):
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

        motifs = greedy_motif_search(dnas, k)
        self.assertTrue(AreStringListsEqual(expected_motifs, motifs))

    def test_extra_dataset(self):
        self._test('Datasets/GreedyMotifSearch/extra.txt')

    def test_sample_dataset(self):
        self._test('Datasets/GreedyMotifSearch/sample.txt')

    def test_dataset_1(self):
        self._test('Datasets/GreedyMotifSearch/data01.txt')

    def test_dataset_2(self):
        self._test('Datasets/GreedyMotifSearch/data02.txt')

    def test_dataset_3(self):
        self._test('Datasets/GreedyMotifSearch/data03.txt')

    def test_dataset_4(self):
        self._test('Datasets/GreedyMotifSearch/data04.txt')


class GreedyMotifSearchLaplace(unittest.TestCase):
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

        motifs = greedy_motif_search_laplace(dnas, k)
        self.assertTrue(AreStringListsEqual(expected_motifs, motifs))

    def test_extra_dataset(self):
        self._test('Datasets/GreedyMotifSearch2/extra.txt')

    def test_sample_dataset(self):
        self._test('Datasets/GreedyMotifSearch2/sample.txt')

    def test_dataset_1(self):
        self._test('Datasets/GreedyMotifSearch2/data01.txt')

    def test_dataset_2(self):
        self._test('Datasets/GreedyMotifSearch2/data02.txt')

    def test_dataset_3(self):
        self._test('Datasets/GreedyMotifSearch2/data03.txt')
