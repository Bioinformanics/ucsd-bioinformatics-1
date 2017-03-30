from random import randint
from Week3.week3_utility import profile_laplace, profile_most_probable_kmer

def score(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]
    max_count = sum([max([c.count(nucleotide) for nucleotide in 'ACGT']) for c in columns])
    return len(motifs[0])*len(motifs) - max_count


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
    rand_ints = [randint(0, len(dna_list[0]) - k) for i in range(t)]
    motifs = [dna_list[i][r:r + k] for i, r in enumerate(rand_ints)]
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
    for repeat in range(1000):
        bms, bm = _randomized_motif_search_cycle(dna_list, k, t)
        if bms < best_motifs_score:
            best_motifs = bm
            best_motifs_score = bms

    return best_motifs


"""
GibbsSampler(Dna, k, t, N)
    randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    BestMotifs ← Motifs
    for j ← 1 to N
        i ← Random(t)
        Profile ← profile matrix constructed from all strings in Motifs except for Motifi
        Motifi ← Profile-randomly generated k-mer in the i-th sequence
        if Score(Motifs) < Score(BestMotifs)
            BestMotifs ← Motifs
    return BestMotifs
"""
def _gibbs_sampler_cycle(dna_list, k, t, N):
    dna_length = len(dna_list[0])
    random_integers = [randint(0, dna_length-k) for i in range(t)]
    motifs = [dna_list[i][r:r+k] for i,r in enumerate(random_integers)]
    best_motifs = list(motifs)
    best_motifs_score = score(best_motifs)
    for j in range(N):
        i = randint(0, t-1)
        profile = profile_laplace([motif for index,motif in enumerate(motifs) if index != i])
        motifs[i] = profile_most_probable_kmer(dna_list[i], k, profile)
        motifs_score = score(motifs)
        if motifs_score < best_motifs_score:
            best_motifs_score = motifs_score
            best_motifs = motifs
    return best_motifs_score, best_motifs

def gibbs_sampler(dna_list, k, t, N):
    best_motifs_score = k*t
    best_motifs = None
    for repeat in range(30):
        bms, bm = _gibbs_sampler_cycle(dna_list, k, t, N)
        if bms < best_motifs_score:
            best_motifs = bm
            best_motifs_score = bms

    return best_motifs

