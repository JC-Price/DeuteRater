# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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


import pymzml  # noqa: 401
import warnings  # noqa: 401
from math import floor

from utils.exc import OutOfMZMLTimeBoundsWarning  # noqa: 401


def retention_time_search(mzml_fp, index_to_ID, search_rt):
    keys = [*index_to_ID]
    left = min(keys)
    right = max(keys)
    mid = 0
    mid_rt = 0

    if (search_rt < mzml_fp[index_to_ID[left]].scan_time_in_minutes() or
            search_rt > mzml_fp[index_to_ID[right]].scan_time_in_minutes()):
        # warnings.warn(OutOfMZMLTimeBoundsWarning(
        #     '{} not within mzml\'s retention time span'.format(
        #         search_rt
        #     )
        # ))
        return -1
    else:
        # TODO: loop problem, check each chunk
        while right - left > 1:
            mid = floor((left + right) / 2)
            # mid_rt = None
            # while mid_rt is None:
            # try:
            mid_rt = mzml_fp[index_to_ID[mid]].scan_time_in_minutes()
            # except Exception:

            if mid_rt < search_rt:
                left = mid
            elif mid_rt > search_rt:
                right = mid
            else:
                return mid
        return right


def get_bounds(mzml_fp, index_to_ID):
    keys = [*index_to_ID]
    mzml_idx_min = min(keys)
    mzml_idx_max = max(keys)
    mzml_rt_min = \
        mzml_fp[index_to_ID[mzml_idx_min]].scan_time_in_minutes()
    mzml_rt_max = \
        mzml_fp[index_to_ID[mzml_idx_max]].scan_time_in_minutes()
    return {
        'idx_min': mzml_idx_min,
        'idx_max': mzml_idx_max,
        'rt_min': mzml_rt_min,
        'rt_max': mzml_rt_max
    }


# def mz_search():
#     pass
