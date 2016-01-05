import pandas as pd

import ccmpred.triplets


def pick_random_triplets(x_pair, n_triplets, min_separation):
    import random
    ncol = x_pair.shape[0]

    out = []
    for _ in range(n_triplets):
        i = random.randint(0, ncol - 2 * min_separation - 1)
        j = random.randint(i + min_separation, ncol - min_separation - 1)
        k = random.randint(j + min_separation, ncol - 1)

        a = random.randint(0, 19)
        b = random.randint(0, 19)
        c = random.randint(0, 19)

        out.append((i, j, k, a, b, c, 0.0))

    return pd.DataFrame(out, columns=("i", "j", "k", "a", "b", "c", "score"), )


# pick 6-triplets ijkabc that maximize w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
def pick_best_ijkabc(x_pair, n_triplets, min_separation):
    return ccmpred.triplets.find_triplet6(x_pair, n_triplets, min_separation)


# pick 3-triplets ijk that maximize sum_{a,b,c=1}^{20} w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
def pick_best_ijk(x_pair, n_triplets, min_separation):
    return ccmpred.triplets.find_triplet3(x_pair, n_triplets, min_separation)


STRATEGIES = {
    'random': pick_random_triplets,
    'best-ijk': pick_best_ijk,
    'best-ijkabc': pick_best_ijkabc,
}

PAIR_TRANSFORMS = {
    'identity': lambda x: x,
    'abs': lambda x: abs(x),
    'square': lambda x: x * x,
}
