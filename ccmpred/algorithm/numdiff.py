import numpy as np
import random


random.seed(42)


def duplicate(x):
    return [np.copy(el) for el in x]


def nd_single(x0):
    ncol = x0[0].shape[0]
    i = random.randint(0, ncol - 1)
    a = random.randint(0, 19)

    return 0, ((i, a), )


def nd_pair(x0):

    if len(x0[1].shape) == 4:

        ncol = x0[0].shape[0]
        i = random.randint(0, ncol - 1)
        j = random.randint(0, ncol - 1)
        a = random.randint(0, 20)
        b = random.randint(0, 20)

        return 1, (
            (i, j, a, b),
            (j, i, b, a)
        )

    elif len(x0[1].shape) == 3:

        ncol = x0[0].shape[0]
        nij = ncol * (ncol - 1) / 2
        ij = random.randint(0, nij - 1)
        a = random.randint(0, 20)
        b = random.randint(0, 20)

        return 1, (
            (ij, a, b),
        )


def nd_triplet(x0):
    ntriplets = x0[2].shape[0]
    t = random.randint(0, ntriplets - 1)

    return 2, ((t, ), )


def numdiff(objfun, x, epsilon=1e-5):
    _, g0 = objfun.evaluate(x)

    x0 = objfun.linear_to_structured(x)
    g0 = objfun.linear_to_structured(g0)

    nd_choices = [nd_single, nd_pair, nd_triplet][:len(x0)]

    print("Pos                                    x                 g            DeltaG")
    while True:
        pos0, positions = random.choice(nd_choices)(x0)

        xA = duplicate(x0)
        xB = duplicate(x0)

        for pos in positions:
            xA[pos0][pos] -= epsilon
            xB[pos0][pos] += epsilon

        fxA, _ = objfun.evaluate(objfun.structured_to_linear(*xA))
        fxB, _ = objfun.evaluate(objfun.structured_to_linear(*xB))

        numdiff = (fxB - fxA) / (2 * epsilon)

        for pos in positions:
            xval = x0[pos0][pos]
            symdiff = g0[pos0][pos]

            posstr = "x{0}[{1}]".format(pos0 + 1, ", ".join("{0}".format(p) for p in pos))
            print("{posstr:20s}   {xval: .10e} {symdiff: .10e}".format(posstr=posstr, xval=xval, symdiff=symdiff,))

        print("gNumeric                                 {numdiff: .10e} {delta: .10e}".format(posstr=posstr, xval=xval, numdiff=numdiff, delta=symdiff - numdiff))

        print()
