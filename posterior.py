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
    parser.add_argument("sum_triplets", help="Triplet sum score file", nargs="?")
    parser.add_argument("out_matrix", help="Output posterior matrix to write")

    parser.add_argument("-m", "--model", default="data/triplet_posterior.pickle", help="mixem probability model learned for computing the posterior")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    sum_matrix = np.loadtxt(args.sum_matrix)

    triplets = collections.defaultdict(list)

    sum_triplets = pd.read_table(args.sum_triplets, header=None, comment="#")
    assert(sum_triplets.shape[1] == 4)
    sum_triplets.columns = ['i', 'j', 'k', 'score']
    for _, trp in sum_triplets.iterrows():
        triplets[(trp.i, trp.j)].append(trp.score)
        triplets[(trp.j, trp.k)].append(trp.score)
        triplets[(trp.i, trp.k)].append(trp.score)

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

            for trp in triplets[i, j]:
                log_bf += (
                    np.log(mixem.probability(trp, *triplet_model[0])) -
                    np.log(mixem.probability(trp, *triplet_model[1]))
                )

            posterior_matrix[i, j] = log_bf
            posterior_matrix[j, i] = log_bf

    posterior_matrix = np.exp(posterior_matrix)
    posterior_matrix = posterior_matrix / (posterior_matrix + 1)

    np.savetxt(args.out_matrix, posterior_matrix)


if __name__ == "__main__":
    main()
