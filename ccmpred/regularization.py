# coding: utf-8
import numpy as np


class L2(object):
    """L2 regularization on single and pair emission potentials"""

    def __init__(self, lambda_single, lambda_pair, center_x_single=0, lambda_triplet=None):
        self.lambda_single = lambda_single
        self.lambda_pair = lambda_pair
        self.center_x_single = center_x_single
        self.lambda_triplet = lambda_triplet

        self.lambda_single_factor = 1
        self.lambda_pair_factor = 1
        self.lambda_triplet_factor = 1

    def __call__(self, x_single, x_pair, x_triplet=None):
        x_ofs = x_single - self.center_x_single

        g_single = 2 * self.lambda_single * x_ofs
        g_pair = 2 * self.lambda_pair * x_pair

        fx_reg = self.lambda_single * self.lambda_single_factor * np.sum(x_ofs * x_ofs) + self.lambda_pair * self.lambda_pair_factor * np.sum(x_pair * x_pair)

        if x_triplet is not None:
            g_triplet = 2 * self.lambda_triplet * x_triplet
            fx_reg += self.lambda_triplet * self.lambda_triplet_factor * np.sum(x_triplet * x_triplet)

            return fx_reg, g_single, g_pair, g_triplet
        else:
            return fx_reg, g_single, g_pair

    def __repr__(self):
        return "L₂ regularization λsingle={0} λpair={1}".format(self.lambda_single, self.lambda_pair)
