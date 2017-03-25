from Week3.week3_utility import *
from utilities import ToSingleLineOfString, ConvertTextToMatrix

def motif_enumeration_problem():
    k = 5
    d = 1
    dnas = ["CGTTACGGAAGTTAAGTATGGGTCG",
            "CGAAAAAGCCTCAGCTTAAACCCAA",
            "GCGTTCTCACCTATACGAAAAGGAA",
            "CGGAACTTCGAAAGACTATAGGTGT",
            "TTTCGTCATTCGTAAGGTGACCTCT",
            "CGCAATGGTGTTCATAAGCGCTTTT"]
    motifs = motif_enumeration(dnas, k, d)

    print(ToSingleLineOfString(motifs))


def calculate_motif_entropy_scores():
    original = ["T   C   G   G   G   G   g   T   T   T   t   t",
                "c   C   G   G   t   G   A   c   T   T   a   C",
                "a   C   G   G   G   G   A   T   T   T   t   C",
                "T   t   G   G   G   G   A   c   T   T   t   t",
                "a   a   G   G   G   G   A   c   T   T   C   C",
                "T   t   G   G   G   G   A   c   T   T   C   C",
                "T   C   G   G   G   G   A   T   T   c   a   t",
                "T   C   G   G   G   G   A   T   T   c   C   t",
                "T   a   G   G   G   G   A   a   c   T   a   C",
                "T   C   G   G   G   t   A   T   a   a   C   C"]
    motifs = map(lambda s: s.upper(), map(lambda s: s.replace(' ', ''), original))
    print(str(sum(score(motifs))))


def medium_string_problem():
    dnas = ["TTCCTGCGTTCCGCACGGCCATGCGTATAAAGCGGTCCTGTC",
            "GCTCGTATCTTTAATTATAACCCGCGGGCGCTCTTCAGCCGT",
            "AATCTGTACATGCCGTTTCGACTTCGGCATCACAAAAGCCGT",
            "GTCGCTACCCCAATGCAAAGCAGTCCTGACAGCAGGCCAGGG",
            "TACGCGTCGCCGAGCTGTGCCGCTAAACCTTATATTCGCACG",
            "AGTCCGACTCCTGCCTGACCGAGCTTATGGAGCTGTGATATC",
            "CTAGGGTGCTTTTACCTTAACCGATGTGAGCCGACGAGCCGT",
            "AGCAGTCAGCTACTAGGTAGCAATCACATACATCGTGGTGAA",
            "AGCGGTCAGCTCCTGATTGCTGCATAAGGAAAGTGATAAGCG",
            "AGCAGTGGCTCACCAGCAAGTAGCTCTGGCGTACCTCGTAAG"]
    k = 6
    ms = find_median_string(dnas, k)
    print(ms)


def profile_most_Probable_k_mer_problem():
    with open('Datasets/ProfileMostProbableKMer/quiz.txt', 'r') as datafile:
        dna = datafile.readline().strip()
        k = int(datafile.readline().strip())
        matrix_text = []
        for loop in range(4):
            matrix_text.append(datafile.readline().strip())
    matrix = ConvertTextToMatrix(matrix_text, float)
    print(get_profile_most_probable_k_mer(dna, matrix))


profile_most_Probable_k_mer_problem()