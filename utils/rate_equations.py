'''Equations

Holds the different equations we can use to fit the data

TODO: Need to know the exact form of the formulas

'''

import numpy as np


def simple(t, k, a, p_adj):
    return a - a * (np.exp(-(k + p_adj) * t))
