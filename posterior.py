#!/usr/bin/env python
import mixem
import pickle
import collections

import numpy as np
import pandas as pd


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("sum_matrix", help="Pair sum score file")
    parser.add_argument("triplets", help="Triplet score file", nargs="?")
    parser.add_argument("out_matrix", help="Output posterior matrix to write")

    parser.add_argument("-m", "--model", default="data/triplet_posterior.pickle", help="mixem probability model learned for computing the posterior")

    args = parser.parse_args()

    return args


def source_triplet(path):
    try:
        mat = pd.read_table(path, sep="\t", comment="#", header=None)
    except:
        return None

    if len(mat.columns) not in (4, 7):
        return None

    elif len(mat.columns) == 4:
        mat.columns = ['i', 'j', 'k', 'confidence']

    elif len(mat.columns) == 7:
        mat.columns = ['i', 'j', 'k', 'a', 'b', 'c', 'confidence']


    return mat

def load_triplets(tripletfile):
    triplets = source_triplet(tripletfile)

    triplet_scores = collections.defaultdict(float)
    for _, trp in triplets.iterrows():
        i, j, k = (int(t) for t in trp[['i', 'j', 'k']])

        triplet_scores[i, j] = max(triplet_scores[i, j], trp.confidence ** 2)
        triplet_scores[j, k] = max(triplet_scores[j, k], trp.confidence ** 2)
        triplet_scores[i, k] = max(triplet_scores[i, k], trp.confidence ** 2)

    return triplet_scores


def main():
    args = parse_args()

    sum_matrix = np.loadtxt(args.sum_matrix)

    triplets = None
    if args.triplets:
        triplets = load_triplets(args.triplets)


    with open(args.model, "rb") as f_model:
        model = pickle.load(f_model)

        pair_model = model['pair_parameters']
        triplet_model = model.get('triplet_parameters', None)

        p_contact = model['p_contact']

    posterior_matrix = np.zeros_like(sum_matrix)

    ncol = sum_matrix.shape[0]

    for i in range(ncol):
        for j in range(i + 1, ncol):

            log_bf = np.log(p_contact) - np.log(1 - p_contact)

            log_bf += (
                np.log(mixem.probability(sum_matrix[i, j], *pair_model[0])) -
                np.log(mixem.probability(sum_matrix[i, j], *pair_model[1]))
            )

            if triplets is not None and (i, j) in triplets:
                log_bf += (
                    np.log(mixem.probability(triplets[(i, j)], *triplet_model[0])) -
                    np.log(mixem.probability(triplets[(i, j)], *triplet_model[1]))
                )

            posterior_matrix[i, j] = log_bf
            posterior_matrix[j, i] = log_bf

    posterior_matrix = np.exp(posterior_matrix)
    posterior_matrix = posterior_matrix / (posterior_matrix + 1)

    np.savetxt(args.out_matrix, posterior_matrix)


if __name__ == "__main__":
    main()
