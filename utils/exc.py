# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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

'''
Module containing all of Deuterater's custom errors and warnings

'''
import warnings  # noqa: F401


class InvalidHeaderError(RuntimeError):
    '''Exception thrown when headers do not match Deuterater's specifications

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidHeaderError, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


class PeakIndexError(IndexError):
    '''Exception thrown when headers do not match Deuterater's specifications

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(PeakIndexError, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


class EmptyIdChunkWarning(RuntimeWarning):
    '''Warning raised when a chunk has no ids

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(EmptyIdChunkWarning, self).__init__(msg + (': %s' % original_e))
        self.original_exception = original_e


class OutOfMZMLTimeBoundsWarning(RuntimeWarning):
    '''Warning raised when retention time is not in the span of the mzml

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(OutOfMZMLTimeBoundsWarning, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e


class InvalidSettingsError(RuntimeError):
    '''Exception thrown when a setting has been entered incorrectly

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidSettingsError, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e


class InvalidSettingsWarning(RuntimeWarning):
    '''Exception thrown when a setting has been entered incorrectly

    Parameters
    ----------
    msg : str
        Human readable string describing the exception.
    original_e : :obj:`Exception`
        original exception thrown, in order to carry through

    '''
    def __init__(self, msg, original_e=None):
        super(InvalidSettingsWarning, self).__init__(
            msg + (': %s' % original_e)
        )
        self.original_exception = original_e
