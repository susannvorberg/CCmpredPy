#!/usr/bin/env python
import argparse
import sys

from ccmpred import AMINO_ACIDS
import ccmpred.io.alignment
import ccmpred.distance
import ccmpred.raw
import ccmpred.counts

import numpy as np

import Bio.PDB

from ccmpred.triplets.pick import STRATEGIES, PAIR_TRANSFORMS


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

    parser.add_argument("-x", "--expand", action="store_true", default=False, help="Expand ijk triplets to all possible abc combinations")
    parser.add_argument("-a", "--amino-letters", action="store_true", default=False, help="Show amino acids as alphabet charactes instead of indices")

    args = parser.parse_args()
    if not args.strategy:
        parser.error("Need to pick a strategy!")

    return args


def main():
    args = parse_args()

    raw = ccmpred.raw.parse(args.rawfile)

    transform = np.vectorize(PAIR_TRANSFORMS[args.pair_transform])

    triplets = STRATEGIES[args.strategy](transform(raw.x_pair), args.n_triplets, args.min_separation)

    if args.expand:
        triplets = ccmpred.triplets.ensure_triplet6(triplets)

    triplets['strategy'] = args.strategy
    triplets['transform'] = args.pair_transform
    triplets['method'] = "pick-{0}-{1}".format(args.strategy, args.pair_transform)

    if args.amino_letters and 'a' in triplets:
        triplets['a'] = triplets['a'].map(lambda a: AMINO_ACIDS[a])
        triplets['b'] = triplets['b'].map(lambda b: AMINO_ACIDS[b])
        triplets['c'] = triplets['c'].map(lambda c: AMINO_ACIDS[c])

    if args.msa_file:
        msa = ccmpred.io.alignment.read_msa(args.msa_file, "psicov")
        assert(msa.shape[1] == raw.ncol)

        if 'a' in triplets.columns:
            tcounts = ccmpred.counts.triplet_counts(msa, ccmpred.triplets.dataframe_to_ndarray(triplets))
            triplets['counts'] = tcounts

        else:
            triplets['counts'] = msa.shape[0]

    if args.pdb_file:
        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure('', args.pdb_file)
        distance_map = ccmpred.distance.distance_map(structure.get_residues())

        assert(distance_map.shape[0] == raw.ncol)
        triplets['dijk'] = distance_map[triplets['i'], triplets['j']] + distance_map[triplets['j'], triplets['k']] + distance_map[triplets['i'], triplets['k']]

    if args.target:
        triplets['target'] = args.target

    triplets.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == '__main__':
    main()
