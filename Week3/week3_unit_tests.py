import unittest

from utilities import ConvertTextToMatrix
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
    def test_sample(self):
        with open('Datasets/ProfileMostProbableKMer/sample.txt', 'r') as datafile:
            dna = datafile.readline().strip()
            k = int(datafile.readline().strip())
            matrix_text = []
            for loop in range(4):
                matrix_text.append(datafile.readline().strip())
        matrix = ConvertTextToMatrix(matrix_text, float)
        k_mer = get_profile_most_probable_k_mer(dna, matrix)
        expected = "CCGAG"
        self.assertEqual(expected, k_mer)

    def test_dataset_1(self):
        with open('Datasets/ProfileMostProbableKMer/data01.txt', 'r') as datafile:
            dna = datafile.readline().strip()
            k = int(datafile.readline().strip())
            matrix_text = []
            for loop in range(4):
                matrix_text.append(datafile.readline().strip())
        matrix = ConvertTextToMatrix(matrix_text, float)
        k_mer = get_profile_most_probable_k_mer(dna, matrix)
        expected = "AGCAGCTT"
        self.assertEqual(expected, k_mer)

    def test_dataset_2(self):
        with open('Datasets/ProfileMostProbableKMer/data02.txt', 'r') as datafile:
            dna = datafile.readline().strip()
            k = int(datafile.readline().strip())
            matrix_text = []
            for loop in range(4):
                matrix_text.append(datafile.readline().strip())
        matrix = ConvertTextToMatrix(matrix_text, float)
        k_mer = get_profile_most_probable_k_mer(dna, matrix)
        expected = "AAGCAGAGTTTA"
        self.assertEqual(expected, k_mer)

    def test_debug_dataset(self):
        with open('Datasets/ProfileMostProbableKMer/debug.txt', 'r') as datafile:
            dna = datafile.readline().strip()
            k = int(datafile.readline().strip())
            matrix_text = []
            for loop in range(4):
                matrix_text.append(datafile.readline().strip())
        matrix = ConvertTextToMatrix(matrix_text, float)
        k_mer = get_profile_most_probable_k_mer(dna, matrix)

        expected = "TGTCGC"
        self.assertEqual(expected, k_mer)



