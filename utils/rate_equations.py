'''Equations

Holds the different equations we can use to fit the data

TODO: Need to know the exact form of the formulas

'''

import numpy as np


def simple(t, k, a, p_adj):
    return a - a * (np.exp(-(k + p_adj) * t))

def calculate_r2 (actual, predicted):
    residuals = actual - predicted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((actual-np.mean(actual))**2)
    r_squared = 1- ss_res/ss_tot
    return r_squared
