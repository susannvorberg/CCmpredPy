import numpy as np
import numpy.ctypeslib as npct
import ctypes
import os.path

array_2d_char = npct.ndpointer(dtype=np.dtype('uint8'), ndim=2, flags='CONTIGUOUS')
array_2d_uint32 = npct.ndpointer(dtype=np.dtype('uint32'), ndim=2, flags='CONTIGUOUS')
array_1d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=1, flags='CONTIGUOUS')
array_2d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=2, flags='CONTIGUOUS')
array_4d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=4, flags='CONTIGUOUS')

libmsac = npct.load_library('libmsacounts', os.path.dirname(__file__))

# void msa_count_single(double *counts, uint8_t *msa, double *weights, uint32_t nrow, uint32_t ncol) {
libmsac.msa_count_single.restype = None
libmsac.msa_count_single.argtypes = [array_2d_float, array_2d_char, array_1d_float, ctypes.c_uint32, ctypes.c_uint32]

# void msa_count_pairs(double *counts, uint8_t *msa, double *weights, uint32_t nrow, uint32_t ncol) {
libmsac.msa_count_pairs.restype = None
libmsac.msa_count_pairs.argtypes = [array_4d_float, array_2d_char, array_1d_float, ctypes.c_uint32, ctypes.c_uint32]

# void msa_count_triplets(double *counts, uint8_t *msa, double *weights, uint32_t nrow, uint32_t ncol, uint32_t ntriplets, uint32_t *triplets) {
libmsac.msa_count_triplets.restype = None
libmsac.msa_count_triplets.argtypes = [array_4d_float, array_2d_char, array_1d_float, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_uint32, array_2d_uint32]

libmsac.msa_char_to_index.restype = None
libmsac.msa_char_to_index.argtypes = [array_2d_char, ctypes.c_uint32, ctypes.c_uint32]

libmsac.msa_index_to_char.restype = None
libmsac.msa_index_to_char.argtypes = [array_2d_char, ctypes.c_uint32, ctypes.c_uint32]


def triplet_counts(msa, triplets, weights=None):

    nrow, ncol = msa.shape
    ntriplets = triplets.shape[0]

    if weights is None:
        weights = np.ones((nrow, ), dtype='float64')

    counts = np.zeros((ntriplets, 21, 21, 21), dtype=np.dtype('float64'))
    libmsac.msa_count_triplets(counts, msa, weights, *msa.shape, ntriplets, triplets)

    return counts


def pair_counts(msa, weights=None):
    nrow, ncol = msa.shape

    if weights is None:
        weights = np.ones((nrow, ), dtype='float64')

    counts = np.zeros((ncol, ncol, 21, 21), dtype=np.dtype('float64'))
    libmsac.msa_count_pairs(counts, msa, weights, *msa.shape)

    return counts


def single_counts(msa, weights=None):
    nrow, ncol = msa.shape

    if weights is None:
        weights = np.ones((nrow, ), dtype='float64')

    counts = np.empty((ncol, 21), dtype=np.dtype('float64'))
    libmsac.msa_count_single(counts, msa, weights, *msa.shape)

    return counts


def both_counts(msa, weights=None):
    pcs = pair_counts(msa, weights)
    scs = np.sum(pcs[np.diag_indices(pcs.shape[0])], axis=1)

    return scs, pcs


def pwm(counts, ignore_gaps=False):
    singles = single_counts(counts)

    if ignore_gaps:
        singles = singles[:, 1:]

    nrow = np.sum(singles, axis=1)

    return singles / nrow[:, np.newaxis]


def index_msa(msa, in_place=False):
    if not in_place:
        msa = msa.copy()

    libmsac.msa_char_to_index(msa, *msa.shape)

    return msa


def char_msa(msa, in_place=False):
    if not in_place:
        msa = msa.copy()

    libmsac.msa_index_to_char(msa, *msa.shape)

    return msa
