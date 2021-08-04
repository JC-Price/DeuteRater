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


# TODO: Eventually look into ctypes or numba for streamlining the memory
#       of this class
# NOTE: This implementation of deuterater does not actually need the 'i'
#       attribute of this class. Consult the program logic in the extract
#       function in extract.py to see why.


class Peak(object):
    '''Contains the data for one peak

    Attributes
    ----------
    mz : float
        The mass-to-charge ratio of the peak
    ab : float
        The abundance value of the peak (also referred to as the intensity)
    i : int
        The peak number (i.e. m0, m1, m2)
    '''
    __slots__ = (
        # defining '__slots__' lets the python interpreter know what fields
        # we will be defining in the class
        'mz',
        'ab',
        'i',
    )

    def __init__(self, mz, abundance, i=None):
        self.mz = mz
        self.ab = abundance
        self.i = i

    # Defining the __repr__ function allows python to call repr()
    # on this object. This is usually much less formatted than the related
    # '__str__' function
    def __repr__(self):
        return 'Peak(mz={0},ab={1},i={2})'.format(
            self.mz, self.ab, self.i
        )

    # Defining the __str__ function allows python to call str()
    # on this object. This is usually the best way to define a 'toString'
    # or similar function
    def __str__(self):
        return 'Peak(mz={:10.2f}, ab={:10.2f}, i={:3d})'.format(
            self.mz, self.ab, self.i
        )
