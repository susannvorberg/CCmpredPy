#!/usr/bin/env python
import argparse

from ccmpred import AMINO_ACIDS
import ccmpred.raw
import ccmpred.triplets


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("rawfile", help="The prediction raw file to use for ranking pairs")
    parser.add_argument("-n", "--n-triplets", type=int, default=1000, help="Maximum number of triplets (default: %(default)d)")
    parser.add_argument("-s", "--min-separation", type=int, default=5, help="Minimum sequence separation (default: %(default)d)")

    args = parser.parse_args()

    return args


def pick_triplets(x_pair, n_triplets=100, min_separation=5):
    triplets = ccmpred.triplets.find_triplets(x_pair, n_triplets, min_separation)

    return triplets


def main():
    args = parse_args()

    raw = ccmpred.raw.parse(args.rawfile)
    triplets = pick_triplets(raw.x_pair, args.n_triplets, args.min_separation)

    for t in range(triplets.shape[0]):
        trp = triplets[t]
        print("[{0:4d}] {1:3d} /{2:3d} /{3:3d}:\t{4}-{5}-{6}".format(t, trp[0], trp[1], trp[2], AMINO_ACIDS[trp[3]], AMINO_ACIDS[trp[4]], AMINO_ACIDS[trp[5]]))

if __name__ == '__main__':
    main()
