#!/usr/bin/env python
import argparse

from ccmpred import AMINO_ACIDS
import ccmpred.io.alignment
import ccmpred.distance
import ccmpred.raw
import ccmpred.triplets
import ccmpred.counts

import numpy as np
import itertools

import Bio.PDB


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

        out.append((i, j, k, a, b, c))

    return np.array(out, dtype="uint32"), np.zeros((n_triplets, ))


# pick 6-triplets ijkabc that maximize w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
def pick_best_ijkabc(x_pair, n_triplets, min_separation):
    return ccmpred.triplets.find_triplet6(x_pair, n_triplets, min_separation)


# pick 3-triplets ijk that maximize sum_{a,b,c=1}^{20} w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
def pick_best_ijk(x_pair, n_triplets, min_separation):
    triplet3, triplet_scores = ccmpred.triplets.find_triplet3(x_pair, n_triplets, min_separation)

    # translate 3-triplets to 6-triplets using (a, b, c) = (A, A, A)
    triplet6 = ccmpred.triplets.triplet3to6(triplet3, short=True)

    triplet_scores = np.repeat(triplet_scores, 20 ** 3)

    return triplet6, triplet_scores


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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("rawfile", help="The prediction raw file to use for ranking pairs")
    parser.add_argument("-m", "--msa-file", default=None, help="An alignment file to count triplets on")
    parser.add_argument("-T", "--pair-transform", default='identity', choices=list(PAIR_TRANSFORMS.keys()), help="Function to transform couplings with")
    parser.add_argument("-t", "--target", default=None, help="Additionally specify a target name to print")
    parser.add_argument("-p", "--pdb-file", default=None, help="Additionally pick a PDB file to determine distance information")
    parser.add_argument("-n", "--n-triplets", type=int, default=1000, help="Maximum number of triplets (default: %(default)d)")
    parser.add_argument("-s", "--min-separation", type=int, default=5, help="Minimum sequence separation (default: %(default)d)")

    parser.add_argument("--strategy", choices=list(STRATEGIES.keys()), default=None, help="Pick a strategy")

    args = parser.parse_args()
    if not args.strategy:
        parser.error("Need to pick a strategy!")

    return args


def main():
    args = parse_args()

    raw = ccmpred.raw.parse(args.rawfile)

    transform = np.vectorize(PAIR_TRANSFORMS[args.pair_transform])

    triplet6, triplet_scores = STRATEGIES[args.strategy](transform(raw.x_pair), args.n_triplets, args.min_separation)

    header_tokens = ['t', 'i', 'j', 'k', 'a', 'b', 'c', 'score']
    header_prefix_tokens = []
    prefix_tokens = []

    header_prefix_tokens.append("strategy")
    prefix_tokens.append(args.strategy)

    header_prefix_tokens.append("transform")
    prefix_tokens.append(args.pair_transform)

    tcounts = None
    if args.msa_file:
        msa = ccmpred.io.alignment.read_msa(args.msa_file, "psicov")
        assert(msa.shape[1] == raw.ncol)
        header_tokens.append("counts")
        tcounts = ccmpred.counts.triplet_counts(msa, triplet6).astype("uint32")

    distance_map = None
    if args.pdb_file:
        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure('', args.pdb_file)
        distance_map = ccmpred.distance.distance_map(structure.get_residues())

        assert(distance_map.shape[0] == raw.ncol)

        header_tokens.append("dijk")

    if args.target:
        header_prefix_tokens.append("target")
        prefix_tokens.append(args.target)

    print("\t".join("{0}".format(tok) for tok in itertools.chain(header_prefix_tokens, header_tokens)))
    for t in range(triplet6.shape[0]):
        trp = triplet6[t]

        tokens = [t, trp[0], trp[1], trp[2], AMINO_ACIDS[trp[3]], AMINO_ACIDS[trp[4]], AMINO_ACIDS[trp[5]], triplet_scores[t]]

        if tcounts is not None:
            tokens.append(tcounts[t])

        if distance_map is not None:
            dist = distance_map[trp[0], trp[1]] + distance_map[trp[0], trp[2]] + distance_map[trp[1], trp[2]]
            tokens.append(dist)

        print("\t".join("{0}".format(tok) for tok in itertools.chain(prefix_tokens, tokens)))

if __name__ == '__main__':
    main()
