#!/usr/bin/env python
import argparse

from ccmpred import AMINO_ACIDS
import ccmpred.io.alignment
import ccmpred.raw
import ccmpred.triplets
import ccmpred.counts


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("rawfile", help="The prediction raw file to use for ranking pairs")
    parser.add_argument("msafile", help="An alignment file to count triplets on")
    parser.add_argument("-n", "--n-triplets", type=int, default=1000, help="Maximum number of triplets (default: %(default)d)")
    parser.add_argument("-s", "--min-separation", type=int, default=5, help="Minimum sequence separation (default: %(default)d)")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    raw = ccmpred.raw.parse(args.rawfile)
    msa = ccmpred.io.alignment.read_msa(args.msafile, "psicov")

    # pick 6-triplets ijkabc that maximize w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
    triplet6, triplet_scores = ccmpred.triplets.find_triplet6(raw.x_pair, args.n_triplets, args.min_separation)
    # # pick 3-triplets ijk that maximize sum_{a,b,c=1}^{20} w_{ij}(a,b) + w_{jk}(b,c) + w_{ik}(a,c)
    # triplet3, triplet_scores = ccmpred.triplets.find_triplet3(raw.x_pair, args.n_triplets, args.min_separation)
    # for t in range(triplet3.shape[0]):
    #     trp = triplet3[t]
    #     print("[{0:4d}] {1:3d} /{2:3d} /{3:3d}".format(t, trp[0], trp[1], trp[2]))

    # # translate 3-triplets to 6-triplets enumerating all abc
    # triplet6 = ccmpred.triplets.triplet3to6(triplet3)

    tcounts = ccmpred.counts.triplet_counts(msa, triplet6)

    print("    t    i    j    k \ta b c  counts        score")
    for t in range(triplet6.shape[0]):
        trp = triplet6[t]
        print("[{0:4d}] {1:3d} /{2:3d} /{3:3d}:\t{4}-{5}-{6}: {7:6.0f} {8:.6e}".format(t, trp[0], trp[1], trp[2], AMINO_ACIDS[trp[3]], AMINO_ACIDS[trp[4]], AMINO_ACIDS[trp[5]], tcounts[t], triplet_scores[t]))

if __name__ == '__main__':
    main()
