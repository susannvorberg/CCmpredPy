#!/usr/bin/env python
import pandas as pd
import Bio.PDB

import ccmpred.distance


def load_triplets(path):

    df = pd.read_table(path)

    if df.shape[1] == 1:

        if len(df.index[0]) == 3:
            df = df.reset_index()
            df.columns = ("i", "j", "k", "score")

        elif len(df.index[0]) == 6:
            df = df.reset_index()
            df.columns = ("i", "j", "k", "a", "b", "c", "score")

    return df[['i', 'j', 'k', 'score']]


def distance_map(pdbfile):
    parser = Bio.PDB.PDBParser()
    resis = list(parser.get_structure('', pdbfile).get_residues())
    return ccmpred.distance.distance_map(resis)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("tripletfile", help="File to load triplets from")
    parser.add_argument("pdbfile", help="File to load PDB file from")
    parser.add_argument("outfile", help="Output file")

    parser.add_argument("-m", "--method", help="Add method information")

    opt = parser.parse_args()

    return opt


def main():
    opt = parse_args()

    triplets = load_triplets(opt.tripletfile)
    distance = distance_map(opt.pdbfile)

    triplets['dijk'] = (
        distance[triplets['i'], triplets['j']] +
        distance[triplets['j'], triplets['k']] +
        distance[triplets['i'], triplets['k']]
    )

    if opt.method:
        triplets['method'] = opt.method

    triplets.to_csv(opt.outfile, sep="\t")


if __name__ == '__main__':
    main()
