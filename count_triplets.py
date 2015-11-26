#!/usr/bin/env python
import numpy as np

import ccmpred.counts
import ccmpred.io.alignment
from ccmpred import AMINO_ACIDS


def main():
    msa = ccmpred.io.alignment.read_msa("data/1atzA.aln", "psicov")

    triplets = np.array([
        (1, 2, 3),
        (2, 3, 4),
        (3, 4, 5)
    ], dtype='uint32')

    triplet_counts = ccmpred.counts.triplet_counts(msa, triplets)

    for t in range(triplet_counts.shape[0]):
        print("\n\nTRIPLET {0}".format(t))

        for ia, a in enumerate(AMINO_ACIDS):

            print("\n{0} @ {1}".format(a, triplets[t, 0]))
            print(" ".join("{0:>4s}".format(aa) for aa in " " + AMINO_ACIDS))

            for ib, b in enumerate(AMINO_ACIDS):
                print("{0:>4s} ".format(b) + " ".join("{0:4.0f}".format(tc) for tc in triplet_counts[t, ia, ib]))

    max_count_pos = np.unravel_index(np.argmax(triplet_counts[:, :20, :20, :20]), (triplet_counts.shape[0], 20, 20, 20))

    print("Highest non-gap count {0} at triplet #{1}: {2} in amino acids {3}".format(triplet_counts[max_count_pos], max_count_pos[0], triplets[max_count_pos[0]], "".join(AMINO_ACIDS[ai] for ai in max_count_pos[1:])))


if __name__ == '__main__':
    main()
