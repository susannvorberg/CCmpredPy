#!/usr/bin/env python
import numpy as np
import argparse
import heapq
import operator


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("matfile", help="The prediction matrix file to use for ranking pairs")

    args = parser.parse_args()

    return args


def pick_triplets(mat, n_triplets=100):
    ncol = mat.shape[0]

    top_triplets = heapq.nlargest(n_triplets, (
        ((i, j, k), mat[i, j] + mat[j, k] + mat[k, i])
        for i in range(ncol)
        for j in range(i + 1, ncol)
        for k in range(j + 1, ncol)
    ), key=operator.itemgetter(1))

    return top_triplets


def main():
    args = parse_args()

    mat = np.loadtxt(args.matfile)
    triplets = pick_triplets(mat)

    print("\n".join(
        "{0}, {1}, {2}: {3}".format(i, j, k, score)
        for (i, j, k), score in triplets
    ))


if __name__ == '__main__':
    main()
