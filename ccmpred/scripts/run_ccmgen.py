#!/usr/bin/env python

import argparse
import os
from ccmpred import CCMpred
import ccmpred.logo
import ccmpred.io.alignment
import ccmpred.raw
import ccmpred.weighting
import ccmpred.sampling
import ccmpred.gaps
import ccmpred.trees
import ccmpred.parameter_handling
import numpy as np

EPILOG = """
Generate a multiple sequence alignment of protein sequences generated from a Markov Random Field model.

In a first step, coupling potentials will have to be learned from a source protein MSA using 
e.g. CCMpredPy with the -b command. 
This can then be passed to the CCMgen call.

"""


def parse_args():
    parser = argparse.ArgumentParser(epilog=EPILOG)

    parser.add_argument("rawfile", help="Raw coupling potential file as generated by the CCMpred --write-msgpack option")
    parser.add_argument("outalnfile", help="Output alignment file for sampled sequences.")


    grp_tr = parser.add_argument_group("Phylogenetic Tree Options")
    grp_tr_me = grp_tr.add_mutually_exclusive_group()
    grp_tr_me.add_argument("--tree-newick",    dest="tree_file", type=str, default=None,
                        help="Load tree from newick-formatted file")
    grp_tr_me.add_argument("--tree-binary",    dest="tree_source", action="store_const", const="binary",
                        help="Generate binary tree")
    grp_tr_me.add_argument("--tree-star",      dest="tree_source", action="store_const", const="star",
                        help="Generate star tree")
    grp_tr_me.add_argument("--mcmc-sampling",      dest="mcmc", action="store_true", default=False,
                        help="Generate MCMC sample without following tree topology.")



    grp_opt = parser.add_argument_group("General Options")
    grp_opt.add_argument("--alnfile", dest="alnfile", type=str, default=None,
                               help="Input alignment file to to specify Neff and NSEQ")
    grp_opt.add_argument("--num-sequences", dest="nseq", metavar="NSEQ", type=int, default=2**10,
                        help="Set the number of sequences to generate to NSEQ "
                             "(does not apply when newick file is specified) [default: 1024]")
    grp_opt.add_argument("--max-gap-pos",  dest="max_gap_pos", default=100, type=int,
                        help="Ignore alignment positions with > MAX_GAP_POS percent gaps. "
                             "[default: %(default)s == no removal of gaps]")
    grp_opt.add_argument("--max-gap-seq",  dest="max_gap_seq",  default=100, type=int,
                        help="Remove sequences with >X percent gaps. "
                             "[default: %(default)s == no removal of sequences]")
    grp_opt.add_argument("--aln-format", dest="aln_format", type=str, default="fasta",
                        help="Specify format for alignment files [default: %(default)s]")
    grp_opt.add_argument("-t", "--num_threads", dest="num_threads", type=int, default=1,
                        help="Specify the number of threads. [default: %(default)s]")



    grp_tr_opt = parser.add_argument_group("Tree Sampling Options")
    grp_tr_opt_me = grp_tr_opt.add_mutually_exclusive_group()
    grp_tr_opt_me.add_argument("--mutation-rate", dest="mutation_rate", type=float, default=1.0,
                        help="Specify constant mutation rate [default: %(default)s]")
    grp_tr_opt_me.add_argument("--mutation-rate-neff", dest="mutation_rate_neff", action="store_true", default=False,
                        help="Set mutation rate to generate alignment with Neff comparable to original MSA "
                             "(requires alignment file with --alnfile).")


    grp_s0 = parser.add_argument_group("Initial Sequence Options")
    grp_s0.add_argument("--seq0-file",  dest="seq0_file",    default=None, type=str, metavar="SEQFILE",
                        help="Specify ancestor sequence in SEQFILE.")
    grp_s0.add_argument("--seq0-mrf",   dest="seq0_mrf", type=int, default=500, metavar="NMUT",
                        help="Sample initial sequence from MRF by mutating a poly-A sequence with NMUT Gibbs steps.")


    grp_mcmc = parser.add_argument_group("MCMC Sampling Options")
    grp_mcmc.add_argument("--mcmc-sample-original", dest="mcmc_sample_type", action="store_const", const="original",
                          default="original", help="Sample sequences starting from original sequences.")
    grp_mcmc.add_argument("--mcmc-sample-random",   dest="mcmc_sample_type", action="store_const", const="random",
                          help="Sample sequences starting from random sequences")
    grp_mcmc.add_argument("--mcmc-sample-random-gapped", dest="mcmc_sample_type", action="store_const",
                          const="random-gapped",
                          help="Sample sequences starting from random sequences but keeping original gap structures "
                               "[default: %(default)s]")
    grp_mcmc.add_argument("--mcmc-burn-in", dest="mcmc_burn_in", type=int, default=500,
                          help="Number of Gibbs sampling steps before a sample is obtained.")




    opt = parser.parse_args()

    if not opt.mcmc and not opt.tree_source and not opt.tree_file:
        parser.error("Need one of the --tree-* options or --mcmc-sampling!")

    if opt.mcmc and not opt.alnfile:
        parser.error("Need alignment file for MCMC sampling!")

    if opt.mutation_rate_neff and not opt.alnfile:
        parser.error("Need alignment file for determining Neff!")

    if opt.tree_source or opt.tree_file:
        if not opt.mutation_rate and not opt.mutation_rate_neff:
            parser.error("Need one of the --mutation-rate* options!")

        if opt.mutation_rate_neff and not opt.seq0_mrf:
            parser.error("Need to specify the number of Gibbs steps for obtaineing seq0 from MRF (--seq0-mrf) when using --mutation_rate_neff !")

    return opt


def main():

    # read command line options
    opt = parse_args()

    ccmpred.logo.logo(what_for="ccmgen")

    # set OMP environment variable for number of threads
    os.environ['OMP_NUM_THREADS'] = str(opt.num_threads)
    print("Using {0} threads for OMP parallelization.".format(os.environ["OMP_NUM_THREADS"]))

    # instantiate CCMpred
    ccm = CCMpred()

    # specify possible file paths
    ccm.set_initraw_file(opt.rawfile)
    ccm.set_pdb_file(opt.pdbfile)


    # read alignment and remove gapped sequences and positions
    if opt.alnfile:
        ccm.set_alignment_file(opt.alnfile)
        ccm.read_alignment(opt.aln_format, opt.max_gap_pos, opt.max_gap_seq)


    #read potentials from binary raw file (possibly remove positions with many gaps)
    ccm.intialise_potentials()
    x = ccmpred.parameter_handling.structured_to_linear(ccm.x_single, ccm.x_pair, nogapstate=True, padding=False)
    ncol = ccm.x_single.shape[0]


    #if MCMC sampling is specified (requires alignment file)
    if opt.mcmc:
        msa_sampled, neff = ccmpred.sampling.generate_mcmc_sample(
            x, ccm.msa, size=opt.nseq, burn_in=opt.mcmc_burn_in, sample_type=opt.mcmc_sample_type)

        ids = ["seq {0}".format(i) for i in range(msa_sampled.shape[0])]

    else:

        tree = ccmpred.trees.CCMTree()

        #prepare tree topology
        if opt.tree_file:
            tree.load_tree(opt.tree_file)
        elif opt.tree_source is not None:

            tree.specify_tree(opt.nseq, opt.tree_source)

        ids = tree.ids

        # sample alignment with Neff similar to alignment Neff (requires alignment file and burn-in)
        if opt.mutation_rate_neff:
            msa_sampled, neff = ccmpred.sampling.sample_to_neff_increasingly(
                tree, ccm.neff_entropy, ncol, x, opt.seq0_mrf)
        # sample alignment with specified mutation rate
        elif opt.mutation_rate:
            if opt.seq0_mrf:
                seq0 = ccmpred.trees.get_seq0_mrf(x, ncol, opt.seq0_mrf)
                print("Ancestor sequence (polyA --> {0} gibbs steps --> seq0) : {1}".format(
                    opt.seq0_mrf, "".join([ccmpred.io.alignment.AMINO_ACIDS[c] for c in seq0[0]])))
            elif opt.seq0_file:
                seq0 = ccmpred.io.alignment.read_msa(opt.seq0_file, opt.aln_format)
                print("Ancestor sequence: {0}".format("".join([ccmpred.io.alignment.AMINO_ACIDS[c] for c in seq0[0]])))
            else:
                seq0 = np.zeros((1, ncol), dtype="uint8")

            msa_sampled, neff = ccmpred.sampling.sample_with_mutation_rate(
                tree, seq0, x, opt.mutation_rate)



    # if gappy positions have been removed
    # insert columns with gaps at that position
    if ccm.max_gap_pos < 100:
        msa_sampled = ccmpred.gaps.backinsert_gapped_positions_aln(
            msa_sampled, ccm.gapped_positions
        )


    print("\nWriting sampled alignment to {0}".format(opt.outalnfile))
    with open(opt.outalnfile, "w") as f:
        descs=["synthetic sequence generated with CCMgen" for _ in range(msa_sampled.shape[0])]
        ccmpred.io.alignment.write_msa(f, msa_sampled, ids, is_indices=True, format=opt.aln_format, descriptions=descs)


if __name__ == '__main__':
    main()
