import operator
import collections

from ccmpred.raw.ccmraw import stream_or_file

_SCORE_COMBINATORS = {}


def score_combinator(name):
    """Decorator to register a score combination function"""
    def inner(f):
        _SCORE_COMBINATORS[name] = f
        return f
    return inner


def get_combinator(name):
    return _SCORE_COMBINATORS[name]


@score_combinator("max")
def combine_max(triplets):
    max_scores = collections.defaultdict(float)
    for (i, j, k, a, b, c), score in triplets:
        max_scores[(i, j, k)] = max(max_scores[(i, j, k)], score)

    return list(max_scores.items())


@score_combinator("square")
def combine_squaresum(triplets):
    sum_scores = collections.defaultdict(float)
    for (i, j, k, a, b, c), score in triplets:
        sum_scores[(i, j, k)] += score ** 2

    return list(sum_scores.items())


@stream_or_file("w")
def write_triplets(f, raw):
    if 'triplets' not in raw.extra_results:
        raise Exception("Raw data does not have triplets information!")

    if 'x_triplet' not in raw.extra_results:
        raise Exception("Raw data does not have triplet potentials!")

    triplets = raw.extra_results['triplets']
    x_triplet = raw.extra_results['x_triplet']

    triplets = list(zip(triplets, x_triplet))

    triplets.sort(key=operator.itemgetter(1), reverse=True)

    f.write("# {0}\n".format(len(triplets)))
    for coords, score in triplets:
        f.write("{0}\t{1:.8e}\n".format("\t".join("{0}".format(te) for te in coords), score))


@stream_or_file("w")
def write_combined_triplets(f, raw, combine=combine_squaresum):
    if 'triplets' not in raw.extra_results:
        raise Exception("Raw data does not have triplets information!")

    if 'x_triplet' not in raw.extra_results:
        raise Exception("Raw data does not have triplet potentials!")

    triplets = raw.extra_results['triplets']
    x_triplet = raw.extra_results['x_triplet']

    combined_triplets = combine(zip(triplets, x_triplet))

    combined_triplets.sort(key=operator.itemgetter(1), reverse=True)

    f.write("# {0}\n".format(len(combined_triplets)))
    for coords, score in combined_triplets:
        f.write("{0}\t{1:.8e}\n".format("\t".join("{0}".format(te) for te in coords), score))
