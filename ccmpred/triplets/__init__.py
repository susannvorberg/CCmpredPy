import numpy as np
import pandas as pd
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

    result = pd.DataFrame(triplets[:ntrip, :], columns=("i", "j", "k", "a", "b", "c"))
    result['score'] = triplet_scores[:ntrip]

    return result


def find_triplet3(x_pair, n_triplets, min_separation):
    triplets = np.zeros((n_triplets, 3), dtype="uint32")
    triplet_scores = np.zeros((n_triplets, ), dtype="double")
    ntrip = libtriplets.find_triplet3(x_pair, triplets, triplet_scores, x_pair.shape[0], n_triplets, min_separation)

    result = pd.DataFrame(triplets[:ntrip, :], columns=("i", "j", "k"))
    result['score'] = triplet_scores[:ntrip]

    return result


def triplet3to6(df):

    df = pd.concat((
        pd.concat((
            df,
            pd.DataFrame({
                "a": np.repeat(a, df.shape[0]),
                "b": np.repeat(b, df.shape[0]),
                "c": np.repeat(c, df.shape[0])
            })
        ), axis=1)
        for a in range(20)
        for b in range(20)
        for c in range(20)
    ))

    return df


def ensure_triplet6(df):
    if 'a' not in df.columns:
        return triplet3to6(df)
    else:
        return df


def dataframe_to_ndarray(df):
    return np.array(df[['i', 'j', 'k', 'a', 'b', 'c']], dtype='uint32')
