from random import randint
from Week3.week3_utility import score, profile_laplace, profile_most_probable_kmer

"""
    RandomizedMotifSearch(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs
"""
def _randomized_motif_search_cycle(dna_list, k, t):
    rand_ints = [randint(0,len(dna_list[0])-k) for i in range(t)]
    motifs = [dna_list[i][r:r+k] for i,r in enumerate(rand_ints)]
    best_motifs = motifs
    best_motifs_score = score(best_motifs)
    while True:
        p = profile_laplace(motifs)
        motifs = [profile_most_probable_kmer(dna, k, p) for dna in dna_list]
        motifs_score = score(motifs)
        if motifs_score < best_motifs_score:
            best_motifs = motifs
            best_motifs_score = motifs_score
        else:
            return best_motifs_score, best_motifs

def randomized_motif_search(dna_list, k, t):
    best_motifs_score = k*t
    best_motifs = None
    for repeat in range(10000):
        bms, bm = _randomized_motif_search_cycle(dna_list, k, t)
        if bms < best_motifs_score:
            best_motifs = bm
            best_motifs_score = bms
        else:
            return best_motifs