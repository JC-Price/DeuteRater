# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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


from abc import ABC, abstractmethod

# TODO: determine which columns are required


class BaseConverter(ABC):
    def __init__(self):
        super().__init__()

    @property
    @abstractmethod
    def converted_data(self):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def load_files(cls, *args):
        return NotImplementedError

    @classmethod
    @abstractmethod
    def convert(cls):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def write(cls, *args):
        raise NotImplementedError

    # NOTE: add functions for additional formats here if we need them

    @staticmethod
    def assert_columns(df):
        # TODO: determine which columns are required
        # TODO: how should we apply this assertion?
        required_columns = []
        assert all([col in list(df) for col in required_columns])
