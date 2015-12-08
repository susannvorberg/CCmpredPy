import numpy as np
import numpy.ctypeslib as npct
import ctypes
import os.path

array_4d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=4, flags='CONTIGUOUS')
array_2d_uint32 = npct.ndpointer(dtype=np.dtype('uint32'), ndim=2, flags='CONTIGUOUS')

libtriplets = npct.load_library('libtriplets', os.path.join(os.path.dirname(__file__), '_build'))

libtriplets.find_triplet6.restype = ctypes.c_uint64
libtriplets.find_triplet6.argtypes = [
    array_4d_float,    # *x_pair
    array_2d_uint32,    # *triplets
    ctypes.c_uint32,    # ncol
    ctypes.c_uint32,    # ntriplets
    ctypes.c_uint32,    # min_separation
]

libtriplets.find_triplet3.restype = ctypes.c_uint64
libtriplets.find_triplet3.argtypes = [
    array_4d_float,    # *x_pair
    array_2d_uint32,    # *triplets
    ctypes.c_uint32,    # ncol
    ctypes.c_uint32,    # ntriplets
    ctypes.c_uint32,    # min_separation
]


def find_triplet6(x_pair, n_triplets, min_separation):
    triplets = np.zeros((n_triplets, 6), dtype="uint32")
    ntrip = libtriplets.find_triplet6(x_pair, triplets, x_pair.shape[0], n_triplets, min_separation)
    return triplets[:ntrip, :]


def find_triplet3(x_pair, n_triplets, min_separation):
    triplets = np.zeros((n_triplets, 3), dtype="uint32")
    ntrip = libtriplets.find_triplet3(x_pair, triplets, x_pair.shape[0], n_triplets, min_separation)
    return triplets[:ntrip, :]
