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

from collections.abc import MutableSequence
from copy import deepcopy

try:
    from obs.peak import Peak  # noqa: 401
    from utils.exc import PeakIndexError
except:
    from DeuteRater.obs.peak import Peak
    from DeuteRater.utils.exc import PeakIndexError

# TODO: docstrings KC
# TODO: is this the best name? KC
# TODO: are properties appropriate? KC
# TODO: read in baseline
# TODO: do i need a len_m function? see how many times it would be used?


class Envelope(MutableSequence):
    '''Contains isotoptic envelope data

    TODO: we need to discuss how technical we want this definition to be

    Attributes
    ----------
    _peaks : :obj:`list` of int
        The list of associated isotopic peaks
    _rt : float
        Retention time of the envelope
    _baseline : float
        The estimated baseline at the locality of the envelope
    _nlb : int
        Number of peaks to look back when identifying negative peaks
    _nla : int
        Number of peaks to look ahead when identifying extra peaks
    _is_valid : bool
        Boolean flag of whether the envelope is valid or not
    '''

    # defining '__slots__' lets the python interpreter know what fields
    # we will be defining in the class
    __slots__ = (
        '_peaks',
        '_rt',
        '_baseline',
        '_nlb',
        '_nla',
        '_is_valid',
    )

    def __init__(self, peaks=[], rt=None, n_lookback=None, n_lookahead=None):
        self._peaks = []
        [self.append_peak(peak) for peak in peaks]
        self._rt = rt
        self._baseline = -1
        self._nlb = n_lookback
        self._nla = n_lookahead
        self._is_valid = True
        self.baseline_calculated = False

    def __deepcopy__(self, memodict={}):
        copy_object = Envelope()
        copy_object._peaks = deepcopy(self._peaks)
        copy_object._rt = self._rt
        copy_object._baseline = self._baseline
        return copy_object

    # Defining the __repr__ function allows python to call repr()
    # on this object. This is usually much less formatted than the related
    # '__str__' function
    def __repr__(self):
        return 'Envelope({})'.format(
            ','.join([repr(peak) for peak in self._peaks])
        )

    # Defining the __str__ function allows python to call str()
    # on this object. This is usually the best way to define a 'toString'
    # or similar function
    def __str__(self):
        return 'Envelope({})'.format(
            ',\n\t'.join([str(peak) for peak in self._peaks])
        )

    # Defining a __getitem__ function enables accessing using index notation
    # please look into the python documentation for additional information
    def __getitem__(self, key):
        if isinstance(key, int):
            if not (0 <= key < len(self)):
                raise IndexError(
                    'index out of bounds: {}'.format(key))
            return self._peaks[key]

        elif isinstance(key, slice):
            if not (0 <= key.start < len(self)):
                raise IndexError(
                    'start index out of bounds: {}'.format(key))
            if not (0 < key.stop <= len(self)):
                raise IndexError(
                    'stop index out of bounds: {}'.format(key))
            return self._peaks[key]

        else:
            raise TypeError(
                'Indexing by "{}" not supported'.format(type(key)))

    # Defining a __setitem__ function enables modification using index notation
    # please look into the python documentation for additional information
    def __setitem__(self, key, value):
        if isinstance(key, int):
            if not (0 <= key < len(self)):
                raise IndexError(
                    'index out of bounds: {}'.format(key))
            if not isinstance(value, Peak):
                raise TypeError(
                    '"{}" is not of type Peak'.format(str(value)))
            self._peaks[key] = value

        elif isinstance(key, slice):
            if not (0 <= key.start < len(self)):
                raise IndexError(
                    'start index out of bounds: {}'.format(key))
            if not (0 < key.stop <= len(self)):
                raise IndexError(
                    'stop index out of bounds: {}'.format(key))

            indices = list(range(len(self)))[key]
            if len(indices) != len(value):
                raise RuntimeError(
                    '{} and {} do not contain the same number of elements'.
                    format(str(value), str(slice)))
            for i in range(len(self))[key]:
                self[i] = value
        else:
            raise TypeError(
                'Indexing by "{}" not supported'.format(type(key)))

    # Defining a __delitem__ function enables deletion using index notation
    # please look into the python documentation for additional information
    def __delitem__(self, key):
        if isinstance(key, int):
            if not (0 <= key < len(self)):
                raise IndexError(
                    'index out of bounds: {}'.format(key))
            del self._peaks[key]
        elif isinstance(key, slice):
            if not (0 <= key.start < len(self)):
                raise IndexError(
                    'start index out of bounds: {}'.format(key))
            elif not (0 < key.stop <= len(self)):
                raise IndexError(
                    'stop index out of bounds: {}'.format(key))
            del self._peaks[key]
        else:
            raise TypeError(
                'Indexing by "{}" not supported'.format(type(key)))

    # Defining a __len__ function enables python to call len() on this object
    def __len__(self):
        return len(self._peaks)

    def insert(self, key, value):
        if not isinstance(value, Peak):
            raise TypeError('"{}" is not of type Peak'.format(value))
        if not isinstance(key, int):
            raise TypeError(
                'Indexing by "{}" not supported'.format(type(key)))
        if not (0 <= key < len(self)):
            raise IndexError(
                'index out of bounds: {}'.format(key))
        self._peaks.insert(key, value)

    def insert_m(self, peak_key, peak):
        if peak_key < 0:
            self._nlb += 1
        self.insert(self, peak_key + self._nlb, peak)

    def get_m(self, peak_key):
        try:
            if isinstance(peak_key, slice):
                return self[slice(
                    peak_key.start + self._nlb,
                    peak_key.stop + self._nlb,
                    peak_key.step
                )]
            else:
                return self[peak_key + self._nlb]
        except IndexError as e:
            raise PeakIndexError(
                'Invalid isotope peak index', original_e=e)

    def has_m(self, peak_key):
        for peak in self._peaks:
            if peak.i == peak_key:
                return True
        return False

    def set_m(self, peak_key, value):
        try:
            if isinstance(peak_key, slice):
                self[slice(
                    peak_key.start + self._nlb,
                    peak_key.stop + self._nlb,
                    peak_key.step
                )] = value
            else:
                self[peak_key + self._nlb] = value
        except IndexError as e:
            raise PeakIndexError(
                'Invalid isotope peak index', original_e=e)

    def del_m(self, peak_key):
        try:
            if isinstance(peak_key, slice):
                s = slice(
                    peak_key.start + self._nlb,
                    peak_key.stop + self._nlb,
                    peak_key.step
                )
                self._nlb -= \
                    sum((i - self._nlb) < 0 for i in range(len(self.peaks)[s]))
                del self[s]
            else:
                self._nlb -= 1
                del self[peak_key + self._nlb]
        except IndexError as e:
            raise PeakIndexError(
                'Invalid isotope peak index', original_e=e)

    def append_peak(self, peak):
        if not isinstance(peak, Peak):
            raise TypeError('"{}" is not of type Peak'.format(peak))
        self._peaks.append(peak)

    def prepend_peak(self, peak):
        if not isinstance(peak, Peak):
            raise TypeError('"{}" is not of type Peak'.format(peak))
        self._nlb += 1
        self.insert(0, peak)

    def normalize(self, method='sum'):
        # should we normalize negative peaks?
        if method == 'sum':
            factor = sum([peak.ab for peak
                         in self.get_peaks()])
        elif method == 'zero':
            factor = self.get_m(0)
        # NOTE: Here is the code for finding the max peak if
        #   We don't want to only use m0
        elif method == 'max':
            factor = max([peak.ab for peak
                         in self.get_peaks()])
        elif method == 'none':
            pass
        else:
            raise RuntimeError(
                '"{}" is not a valid normalization method'.format(method))
        if method != 'none':
            for peak in self._peaks:
                peak.ab /= factor
        # self._norms.put(factor)

    def unnormalize(self, n_undo=1):
        for i in range(n_undo):
            for peak in self._peaks:
                peak.ab *= self._norms.get()

    @property
    def is_valid(self):
        return self._is_valid

    @is_valid.setter
    def is_valid(self, val):
        self._is_valid = val

    @property
    def baseline(self):
        return self._baseline

    @baseline.getter
    def baseline(self):
        if type(self._baseline) == list:    # if the baseline is still represeted as a list, then calc the actual baseline.
            from numpy import median
            def mad(values):
                m = median(values)
                return median([abs(a-m) for a in values])
            normal_distribution_scale_factor = 1.4826
            self._baseline = mad(self._baseline) * 3 * normal_distribution_scale_factor
        return self._baseline

    @baseline.setter
    def baseline(self, val):
        self._baseline = val

    @property
    def rt(self):
        return self._rt

    @rt.setter
    def rt(self, val):
        self._rt = val

    def get_peaks(self):
        return self._peaks[self._nlb:-self._nla]

    def get_negative_peaks(self):
        return self._peaks[:self._nlb]

    def get_lookahead_peaks(self):
        return self._peaks[-self._nla:]

    # TODO: Figure out how to make these _obs functions more self-documenting
    # TODO: Perform this in only one pass

    # The obs functions are for more streamlined access to the mz and abundance
    #   values of the peaks in this envelope
    def to_obs(self):
        '''
        Returns
        -------
        :obj:`tuple` of :obj:`tuple` of float
        '''
        return (tuple(peak.mz for peak in self.get_peaks()),
                tuple(peak.ab for peak in self.get_peaks()))

    def lb_obs(self):
        return (tuple(peak.mz for peak in self.get_negative_peaks()),
                tuple(peak.ab for peak in self.get_negative_peaks()))

    def la_obs(self):

        return (tuple(peak.mz for peak in self.get_lookahead_peaks()),
                tuple(peak.ab for peak in self.get_lookahead_peaks()))
