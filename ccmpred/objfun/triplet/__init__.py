import numpy as np

import ccmpred.raw
import ccmpred.regularization
import ccmpred.objfun
import ccmpred.objfun.triplet.cext
import ccmpred.triplets


class TripletPseudoLikelihood(ccmpred.objfun.ObjectiveFunction):

    def __init__(self, msa, freqs, weights, regularization, triplets):
        super(TripletPseudoLikelihood, self).__init__()

        self.msa = msa
        self.triplets = triplets
        self.nrow, self.ncol = msa.shape
        self.weights = weights
        self.regularization = regularization

        # TODO make this configurable
        self.regularization.lambda_triplet = 0.1

        neff = np.sum(weights)
        freqs_single, freqs_pair = freqs

        msa_counts_single, msa_counts_pair = neff * freqs_single, neff * freqs_pair
        msa_counts_triplets = ccmpred.counts.triplet_counts(msa, triplets, weights)

        msa_counts_single[:, 20] = 0
        msa_counts_pair[:, :, 20, :] = 0
        msa_counts_pair[:, :, :, 20] = 0

        for i in range(self.ncol):
            msa_counts_pair[i, i, :, :] = 0

        self.g_init = structured_to_linear(msa_counts_single, i_j_to_ij(msa_counts_pair * 2), msa_counts_triplets * 3)

        self.nvar_single = self.ncol * 21
        self.nvar_pair = self.ncol * (self.ncol - 1) / 2 * 21 * 21
        self.nvar_triplet = triplets.shape[0]
        self.nvar = self.nvar_single + self.nvar_pair + self.nvar_triplet

        # memory allocation for intermediate variables
        self.g = np.empty((self.nvar,), dtype=np.dtype('float64'))

        self.linear_to_structured = lambda x: linear_to_structured(x, self.ncol, self.triplets.shape[0])
        self.structured_to_linear = structured_to_linear

    @classmethod
    def init_from_default(cls, msa, freqs, weights, regularization, strategy, transform):
        res = cls(msa, freqs, weights, regularization)

        # TODO
        raise Exception("TODO - run regular PLL then init this!")

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
    def init_from_raw(cls, msa, freqs, weights, raw, regularization, strategy, transform):

        if msa.shape[1] != raw.ncol:
            raise Exception('Mismatching number of columns: MSA {0}, raw {1}'.format(msa.shape[1], raw.ncol))

        # TODO make n_triplets and min_separation configurable
        n_triplets = 1000
        min_separation = 5

        print("Picking up to {0} triplets with separation {1} using strategy {2} and transform {3}".format(n_triplets, min_separation, strategy, transform))
        triplets = strategy(transform(raw.x_pair), n_triplets, min_separation)

        # no matter what stategy, ensure that triplets is a shape (T, 6) ndarray of triplets
        triplets = ccmpred.triplets.ensure_triplet6(triplets)
        triplets = ccmpred.triplets.dataframe_to_ndarray(triplets)

        res = cls(msa, freqs, weights, regularization, triplets)

        x = structured_to_linear(raw.x_single, i_j_to_ij(raw.x_pair), np.zeros((triplets.shape[0],)))

        return x, res

    def finalize(self, x):
        ncol = self.ncol

        x_single, x_pair, x_triplet = linear_to_structured(x, ncol, self.triplets.shape[0])

        ox_pair = np.zeros((ncol, ncol, 21, 21))
        ij = 0
        for i in range(ncol):
            for j in range(i + 1, ncol):
                ox_pair[i, j] = x_pair[ij]
                ox_pair[j, i] = x_pair[ij]
                ij += 1

        extra_results = {
            "x_triplet": x_triplet,
            "triplets": self.triplets
        }

        return ccmpred.raw.CCMRaw(ncol, x_single[:, :20], ox_pair, {}, extra_results=extra_results)

    def evaluate(self, x):
        fx, g = ccmpred.objfun.triplet.cext.evaluate(x, self.g, self.weights, self.msa, self.triplets)

        self.g -= self.g_init

        x_single, x_pair, x_triplet = linear_to_structured(x, self.ncol, self.triplets.shape[0])
        g_single, g_pair, g_triplet = linear_to_structured(g, self.ncol, self.triplets.shape[0])

        fx_reg, g_single_reg, g_pair_reg, g_triplet_reg = self.regularization(x_single[:, :20], x_pair, x_triplet)

        g_reg = structured_to_linear(g_single_reg, g_pair_reg, g_triplet_reg)
        fx += fx_reg
        g += g_reg

        return fx, g

    def __repr__(self):
        return "Triplet-PLL ({0})".format(self.regularization)


def linear_to_structured(x, ncol, ntriplets):
    """Convert linear vector of variables into multidimensional arrays.

    in linear memory, memory order is x1[i, a], x2[ij, a, i] and x3[t] (dimensions Lx21 + (L*(L-1)/2)x21x21 + T)
    output will have  memory order of x1[i, a], x2[ij, a, b] and x3[t] (dimensions Lx21 + (L*(L-1)/2)x21x21 + T)
    """

    nsingle = ncol * 21
    npair = ncol * (ncol - 1) / 2 * 21 * 21

    x_single = x[:nsingle].reshape((ncol, 21))
    x_pair = x[nsingle:(nsingle + npair)].reshape((ncol * (ncol - 1) / 2, 21, 21))
    x_triplet = x[(nsingle + npair):].reshape((ntriplets, ))

    return x_single, x_pair, x_triplet


def structured_to_linear(x_single, x_pair, x_triplet):
    """Convert structured variables into linear array

    with input arrays of memory order x1[i, a], x2[ij, a, i] and x3[t] (dimensions Lx21 + (L*(L-1)/2)x21x21 + T)
    output will have  memory order of x1[i, a], x2[ij, a, b] and x3[t] (dimensions Lx21 + (L*(L-1)/2)x21x21 + T)
    """

    ncol = x_single.shape[0]

    nsingle = ncol * 21
    npair = ncol * (ncol - 1) / 2 * 21 * 21
    ntriplet = x_triplet.shape[0]
    nvar = nsingle + npair + ntriplet

    x = np.zeros((nvar, ), dtype='float64')

    xo_single = x[:nsingle].reshape(ncol, 21)
    xo_single[:, :20] = x_single[:, :20]

    x[nsingle:(nsingle + npair)] = x_pair.reshape(-1)
    x[(nsingle + npair):] = x_triplet.reshape(-1)

    return x


def i_j_to_ij(x_pair):
    ncol = x_pair.shape[0]
    xo_pair = np.zeros((ncol * (ncol - 1) / 2, 21, 21))
    ij = 0
    for i in range(ncol):
        for j in range(i + 1, ncol):
            xo_pair[ij] = x_pair[i, j]
            ij += 1

    return xo_pair
