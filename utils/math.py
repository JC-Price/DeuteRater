# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2024  Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
All rights reserved.
Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:
    * Redistributions of source code must retain the
      above copyright notice, this list of conditions
      and the following disclaimer.
    * Redistributions in binary form must reproduce
      the above copyright notice, this list of conditions
      and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors
      may be used to endorse or promote products derived
      from this software without specific prior written
      permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


'''Mathematical Utilities

This module contains mathematical functions applied in different parts of
the extractor. These are fairly standard formulas that are not subject to change


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
    previous_err = np.seterr(divide="ignore", invalid="ignore")
    value = vector/np.linalg.norm(vector)
    np.seterr(**previous_err)
    return value


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
