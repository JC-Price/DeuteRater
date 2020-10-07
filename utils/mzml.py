

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
