import numpy as np
import numpy.ctypeslib as npct
import ctypes
import os.path

array_1d_float = npct.ndpointer(dtype=np.dtype('float64'), ndim=1, flags='CONTIGUOUS')
array_1d_uint32 = npct.ndpointer(dtype=np.dtype('uint32'), ndim=1, flags='CONTIGUOUS')
array_2d_char = npct.ndpointer(dtype=np.dtype('uint8'), ndim=2, flags='CONTIGUOUS')

libtriplet = npct.load_library('libtriplet', os.path.join(os.path.dirname(__file__), '_build'))

libtriplet.evaluate_triplet_pll.restype = ctypes.c_double
libtriplet.evaluate_triplet_pll.argtypes = [
    array_1d_float,    # *x
    array_1d_float,    # *g
    array_1d_float,    # *weights
    array_2d_char,      # *msa
    array_1d_uint32,    # *triplets
    ctypes.c_uint32,    # nrow
    ctypes.c_uint32,    # ncol
    ctypes.c_uint32,    # ntriplets
]


def evaluate(x, g, weights, msa, triplets):
    nrow, ncol = msa.shape
    fx = libtriplet.evaluate_triplet_pll(x, g, weights, msa, triplets.reshape(-1), nrow, ncol, triplets.shape[0])
    return fx, g
