import numpy as np

import ccmpred.raw
import ccmpred.regularization
import ccmpred.objfun
import ccmpred.objfun.pll.cext


class TripletPseudoLikelihood(ccmpred.objfun.ObjectiveFunction):

    def __init__(self, msa, freqs, weights, regularization, triplets):
        super(TripletPseudoLikelihood, self).__init__()

        self.msa = msa
        self.triplets = triplets
        self.nrow, self.ncol = msa.shape
        self.weights = weights
        self.regularization = regularization

        neff = np.sum(weights)
        freqs_single, freqs_pair = freqs
        msa_counts_single, msa_counts_pair = neff * freqs_single, neff * freqs_pair
        msa_counts_triplets = ccmpred.counts.triplet_counts(msa, triplets, weights)

        msa_counts_single[:, 20] = 0
        msa_counts_pair[:, :, 20, :] = 0
        msa_counts_pair[:, :, :, 20] = 0
        msa_counts_triplets[:, :, :, 20] = 0
        msa_counts_triplets[:, :, 20, :] = 0
        msa_counts_triplets[:, 20, :, :] = 0

        for i in range(self.ncol):
            msa_counts_pair[i, i, :, :] = 0

        self.g_init = structured_to_linear(msa_counts_single, msa_counts_pair, msa_counts_triplets)

        self.nvar_single = self.ncol * 21
        self.nvar_pair = self.ncol * (self.ncol - 1) / 2
        self.nvar_triplet = triplets.shape[0] * 21 * 21 * 21
        self.nvar = self.nvar_single + self.nvar_pair + self.nvar_triplet

        # memory allocation for intermediate variables
        self.g = np.empty((self.nvar,), dtype=np.dtype('float64'))

        self.linear_to_structured = lambda x: linear_to_structured(x, self.ncol, self.triplets.shape[0])
        self.structured_to_linear = structured_to_linear

    @classmethod
    def init_from_default(cls, msa, freqs, weights, regularization):
        res = cls(msa, freqs, weights, regularization)

        # TODO centering
        if hasattr(regularization, "center_x_single"):
            ncol = msa.shape[1]
            x_pair = np.zeros((ncol, ncol, 21, 21), dtype="float64")

            # TODO what to do about triplets?
            x = structured_to_linear(regularization.center_x_single, x_pair)

        else:
            x = np.zeros((res.nvar, ), dtype=np.dtype('float64'))

        return x, res

    @classmethod
    def init_from_raw(cls, msa, freqs, weights, raw, regularization):
        res = cls(msa, freqs, weights, regularization)

        if msa.shape[1] != raw.ncol:
            raise Exception('Mismatching number of columns: MSA {0}, raw {1}'.format(msa.shape[1], raw.ncol))

        # TODO what to do with triplets?
        x = structured_to_linear(raw.x_single, raw.x_pair)

        return x, res

    def finalize(self, x):
        x_single, x_pair, x_triplet = linear_to_structured(x, self.ncol, self.triplets.shape[0])
        return ccmpred.raw.CCMRaw(self.ncol, x_single[:, :20], x_pair, {}, x_triplet=x_triplet)

    def evaluate(self, x):

        # TODO update evaluate call
        fx, g = ccmpred.objfun.pll.cext.evaluate(x, self.g, self.g2, self.weights, self.msa)

        # TODO check precomputed counts signs
        self.g -= self.g_init

        x_single, x_pair, x_triplet = linear_to_structured(x, self.ncol, self.triplets.shape[0])
        g_single, g_pair, x_triplet = linear_to_structured(g, self.ncol, self.triplets.shape[0])

        fx_reg, g_single_reg, g_pair_reg, g_triplet_reg = self.regularization(x_single, x_pair, x_triplet)

        g_reg = structured_to_linear(g_single_reg, g_pair_reg, g_triplet_reg)
        fx += fx_reg
        g += g_reg

        return fx, g

    def __repr__(self):
        return "Triplet-PLL ({0})".format(self.regularization)


def linear_to_structured(x, ncol, ntriplets):
    """Convert linear vector of variables into multidimensional arrays.

    in linear memory, memory order is x1[i, a], x2[ij, a, i] and x3[t, a, b, c] (dimensions Lx21 + (L*(L-1)/2)x21x21 + Tx21x21x21)
    output will have  memory order of x1[i, a], x2[ij, a, b] and x3[t, a, b, c] (dimensions Lx21 + (L*(L-1)/2)x21x21 + Tx21x21x21)
    """

    nsingle = ncol * 21
    npair = ncol * (ncol - 1) / 2 * 21 * 21

    x_single = x[:nsingle].reshape((21, ncol))
    x_pair = x[nsingle:(nsingle + npair)].reshape((ncol * ncol - 1, 21, 21))
    x_triplet = x[(nsingle + npair):].reshape((ntriplets, 21, 21, 21))

    return x_single, x_pair, x_triplet


def structured_to_linear(x_single, x_pair, x_triplet):
    """Convert structured variables into linear array

    with input arrays of memory order x1[i, a], x2[ij, a, i] and x3[t, a, b, c] (dimensions Lx21 + (L*(L-1)/2)x21x21 + Tx21x21x21)
    output will have  memory order of x1[i, a], x2[ij, a, b] and x3[t, a, b, c] (dimensions Lx21 + (L*(L-1)/2)x21x21 + Tx21x21x21)
    """

    ncol = x_single.shape[0]
    nt = x_triplet.shape[0]

    nsingle = ncol * 21
    npair = ncol * (ncol - 1) / 2 * 21 * 21
    ntriplet = nt * 21 * 21 * 21
    nvar = nsingle + npair + ntriplet

    x = np.zeros((nvar, ), dtype='float64')
    x[:nsingle] = x_single.T.reshape(-1)
    x[nsingle:(nsingle + npair)] = x_pair.reshape(-1)
    x[npair:] = x_triplet.reshape(-1)

    return x
