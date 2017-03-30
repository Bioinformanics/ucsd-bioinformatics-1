from random import randint


def read_dna_list_from_fasta():
    with open("upstream250.txt", "r") as datafile:
        lines = datafile.readlines()
        return [line.strip() for line in lines if not line.startswith(">")]


def profile_laplace(motifs):
    columns = [''.join(dna) for dna in zip(*motifs)]
    return [[float(column.count(base) + 1) / float(len(column) + 4) for base in 'ACGT'] for column in columns]


def score(motifs):
    columns = [''.join(dna) for dna in zip(*motifs)]
    max_count = sum([max([column.count(base) for base in 'ACGT']) for column in columns])
    return len(motifs[0])*len(motifs) - max_count


def get_profile_most_probable_k_mer(dna, k, profile):

    max_probability = -1
    max_k_mer = None
    for index in range(len(dna) - k + 1):
        k_mer = dna[index: index + k]
        probability = 1
        for pos, base in enumerate(k_mer):
            probability *= profile[pos]['ACGT'.index(base)]
        if probability > max_probability:
            max_probability = probability
            max_k_mer = k_mer
    return max_k_mer


def gibbs_sampler_cycle(dna_list, k, N):
    t = len(dna_list)
    dna_length = len(dna_list[0])
    random_integers = [randint(0, dna_length-k) for i in range(t)]
    motifs = [dna_list[i][r:r+k] for i,r in enumerate(random_integers)]
    best_motifs = list(motifs)
    best_motifs_score = score(best_motifs)
    for j in range(N):
        i = randint(0, t-1)
        profile = profile_laplace([motif for index,motif in enumerate(motifs) if index != i])
        motifs[i] = get_profile_most_probable_k_mer(dna_list[i], k, profile)
        motifs_score = score(motifs)
        if motifs_score < best_motifs_score:
            best_motifs_score = motifs_score
            best_motifs = motifs
    return best_motifs_score, best_motifs


def gibbs_sampler(dna_list, k, N):
    best_motifs_score = k*len(dna_list)
    best_motifs = None
    for repeat in range(20):
        bms, bm = gibbs_sampler_cycle(dna_list, k, N)
        if bms < best_motifs_score:
            best_motifs = bm
            best_motifs_score = bms

    return best_motifs


def week5():
    dna_list = read_dna_list_from_fasta()
    k = 20
    N = 1000

    motifs = gibbs_sampler(dna_list, k, N)

