'''Mathematical Utilities

This module contains mathematical functions applied in different parts of
Deuterater. These are fairly standard formulas that are not subject to change

TODO: Discuss which formulas are necessary here and validate that the
formulas supplied are properly implemented.
NOTE: This is probably a good point to start writing the pytest framework

'''
import numpy as np


def inclusive_range(start, stop, step=1):
    return range(start, stop+1, step)


def inclusive_slice(start, stop, step=1):
    return slice(start, stop+1, step)


def find_nearest_index(arr, val):
    idx = np.searchsorted(arr, val)
    if idx == 0:
        return idx
    elif idx == len(arr):
        return idx - 1
    else:
        return min([idx - 1, idx], key=lambda x: abs(arr[x] - val))


def unit_vector(vector):
    '''Calculates the unit vector of a given vector

    Parameters
    ----------
    vector : :obj:`list` of :obj:`float`

    Returns
    ----------
    vector : :obj:`list` of :obj:`float`

    '''
    return vector/np.linalg.norm(vector)


def angle_between(v1, v2):
    '''Calculates the unit vector of a given vector

    Parameters
    ----------
    v1 : :obj:`list` of :obj:`float`
        The first vector to compare
    v2 : :obj:`list` of :obj:`float`
        The second vector to compare

    Returns
    ----------
    :obj:`float`
        Angle between the two vectors, given in radians

    '''
    uv1 = unit_vector(v1)
    uv2 = unit_vector(v2)
    angle = np.arccos(np.dot(uv1, uv2))
    if np.isnan(angle):
        if(uv1 == uv2).all():
            return 0.0
        else:
            return np.pi
    return angle
