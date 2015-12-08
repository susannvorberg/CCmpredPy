import numpy as np
import numpy.ctypeslib as npct
import ctypes
import itertools
import os.path

array_1d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=1, flags='CONTIGUOUS')
array_4d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=4, flags='CONTIGUOUS')
array_2d_uint32 = npct.ndpointer(dtype=np.dtype('uint32'), ndim=2, flags='CONTIGUOUS')

libtriplets = npct.load_library('libtriplets', os.path.join(os.path.dirname(__file__), '_build'))

libtriplets.find_triplet6.restype = ctypes.c_uint64
libtriplets.find_triplet6.argtypes = [
    array_4d_float,    # *x_pair
    array_2d_uint32,    # *triplets
    array_1d_float,     # *triplet_scores
    ctypes.c_uint32,    # ncol
    ctypes.c_uint32,    # ntriplets
    ctypes.c_uint32,    # min_separation
]

libtriplets.find_triplet3.restype = ctypes.c_uint64
libtriplets.find_triplet3.argtypes = [
    array_4d_float,    # *x_pair
    array_2d_uint32,    # *triplets
    array_1d_float,     # *triplet_scores
    ctypes.c_uint32,    # ncol
    ctypes.c_uint32,    # ntriplets
    ctypes.c_uint32,    # min_separation
]


def find_triplet6(x_pair, n_triplets, min_separation):
    triplets = np.zeros((n_triplets, 6), dtype="uint32")
    triplet_scores = np.zeros((n_triplets, ), dtype="double")
    ntrip = libtriplets.find_triplet6(x_pair, triplets, triplet_scores, x_pair.shape[0], n_triplets, min_separation)
    return triplets[:ntrip, :], triplet_scores[:ntrip]


def find_triplet3(x_pair, n_triplets, min_separation):
    triplets = np.zeros((n_triplets, 3), dtype="uint32")
    triplet_scores = np.zeros((n_triplets, ), dtype="double")
    ntrip = libtriplets.find_triplet3(x_pair, triplets, triplet_scores, x_pair.shape[0], n_triplets, min_separation)
    return triplets[:ntrip, :], triplet_scores[:ntrip]


def triplet3to6(triplet3):
    amino_grid = np.array(list(itertools.product(* ([range(20)] * 3))))

    triplet6 = np.empty((triplet3.shape[0] * amino_grid.shape[0], 6), dtype="uint32")
    triplet6[:, :3] = np.repeat(triplet3, amino_grid.shape[0], axis=0)
    triplet6[:, 3:6] = np.repeat(amino_grid, triplet3.shape[0], axis=0)

    return triplet6
