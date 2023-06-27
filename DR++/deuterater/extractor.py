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

import pandas as pd
import numpy as np  # noqa: 401
import pymzml
from tqdm import tqdm  # noqa: 401

from more_itertools import one
from functools import partial
import multiprocessing as mp
from pathlib import Path  # noqa: 401
import os
import traceback
import warnings  # noqa: 401

from math import ceil
from operator import mul

import deuterater.settings as settings
import utils.mzml as dml
import utils.extract as due
from utils.exc import InvalidHeaderError  # noqa: 401


# TODO: Code for Validating headers
# TODO: n_isos to look for should be supplied in the ID file
#   NOTE: this includes consideration of peaks_included_under_mass_cutoff
# TODO: chunk on lines not number of partitions
# TODO: we achieve best efficiency if the number of chunks is a multiple
#       of the number of processors. Does this need to be implemented?

# TODO: if we trim the IDs, what to do with obs trimmed off?
# TODO: change the logic to be based off of pandas chunking?
# TODO: let the extractor handle gzipped mzml files
# TODO: Speedup with cython
# TODO: add check to see if file is writeable
# TODO: only compose dictionary once instead of twice?
# TODO: How to deal with killing only one instance of extractor
#       in a batch processing scenario?
# TODO: pause execution to close file, don't terminate

# NOTE: not all mzmls have contiguous native IDs


class Extractor:  # TODO name change
    '''Class that handles extracting the relevant data from an mzml file.

    TODO: Add more info on what this does/how it differs from extract

    i.e.

    Attributes
    ----------
    n_processors : int
        The number of partitions in which to process the data.
    max_chunk_size : int
        The upper limit of how large a chunk can be
    ids : :obj:`pandas.Dataframe`
        A dataframe containing identifications made by database software
        on unlabeled mass spectrometry data
    mzml : :obj:`pymzml.run.Reader`
        An object providing access to the mass spec data to be searched
    model : :obj:`pandas.Dataframe`
        The aggregation of the data extracted from `ids` and `mzml`

    '''

    def __init__(self, settings_path, id_path, mzml_path, out_path):
        '''
        Parameters
        ----------
        ids : str
            The name of the file containing the identifications. This data
            will likely have been taken from an unlabele-
            For addition information, see the settings file's documentation

        '''
        self.settings_path = settings_path
        settings.load(self.settings_path)

        self.id_path = id_path
        self.mzml_path = mzml_path
        self.out_path = out_path

        self.ids = {}
        self._id_chunks = ()
        self._mzml_native_id_bounds = []
        self.model = pd.DataFrame()

        try:
            if settings.recognize_available_cores is True:
                # BD: Issue with mp.cpu_count() finding too many cores available
                self._n_processors = round(mp.cpu_count() * 0.75)
                # self._n_processors = mp.cpu_count()
            else:
                self._n_processors = settings.n_processors
            # $breaks windows/python interactions if too many cores are used. very niche application but still relevant
            if self._n_processors > 60:  # changed from 60 back to 60, -Coleman
                self.n_processors = 60
            self._chunk_size = settings.chunk_size
            self._chunking_threshold = mul(
                self._n_processors,
                settings.chunk_size
            )
            self._id_rt_unit = settings.id_file_rt_unit
            self._trim_ids = settings.trim_ids_to_mzml_bounds
            self._rt_window = settings.time_window
            self._mp_pool = mp.Pool(self._n_processors)

            if not os.path.exists(self.out_path):
                open(self.out_path, 'w').close()
            # TODO: this needs fixed, it doesn't catch open files
            if not os.access(self.out_path, os.W_OK):
                raise PermissionError('Output path not writeable.')

        except Exception as e:
            print(e)
            traceback.print_tb(e.__traceback__)
            raise

    def load(self):
        """Loads the associated files into memory

        This is called separately from the initialization of the Extractor
        class so that if many files are given, we only have one Extractor
        fully loaded into memory at a time
        """
        try:
            self._load_ids(self.id_path)
            self._get_mzml_bounds(self.mzml_path)
            self._partition_ids(trim=self._trim_ids)
        except Exception as e:
            print(e)
            traceback.print_tb(e.__traceback__)
            raise

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    def _load_ids(self, filename):
        '''Loads the id file into memory

        Parameters
        ----------
        filename : str
            The name of the file containing the identifications

        Returns
        -------
        :obj:`bool`
            True if the file loads successfullly, False otherwise.

        Raises
        ------
        FileNotFoundError
            If the file does not exist
        InvalidHeaderError
            If the file does not contain the appropriate information.

        '''
        self.ids = pd.read_csv(self.id_path)
        if self._id_rt_unit == 'sec':
            self.ids['rt'] = self.ids['Precursor Retention Time (sec)'].apply(lambda x: x / 60.0)
        elif self._id_rt_unit == 'min':
            self.ids['rt'] = self.ids['Precursor Retention Time (sec)']
        self.ids.sort_values(by=['rt'], inplace=True)
        self.ids.reset_index(inplace=True, drop=True)

        # self.ids['n_isos'] = 5  # TODO: Temp value, see note at top of files
        def num_peaks_by_mass(mass):
            if mass < 1500:
                return 3
            elif mass < 2400:
                return 4
            else:
                return 5

        self.ids['n_isos'] = self.ids['neutromers_to_extract'].astype(
            np.int8)  # TODO: Temp value, see note at top of files
        # self.ids['n_isos'] = self.ids['Peptide Theoretical Mass'].apply(num_peaks_by_mass)

    def _partition_ids(self, trim):
        '''splits the id file

        Returns
        -------
        :obj:`list` of :obj:`pandas.Dataframe`

        '''
        # NOTE: there is not really a reason to chunk the read operation
        # until the ID files are gigabytes in size
        if trim:
            mask = (self._mzml_rt_min - self._rt_window < self.ids['rt']) \
                   & (self.ids['rt'] + self._rt_window < self._mzml_rt_max)
            self.ids = self.ids[mask]
            self.ids.reset_index(inplace=True, drop=True)

        if len(self.ids.index) <= self._chunking_threshold:
            self._chunk_size = ceil(len(self.ids.index) / self._n_processors)

        try:
            self._id_chunks = [
                self.ids.loc[i:i + self._chunk_size - 1, :].copy() for i in range(
                    0, len(self.ids.index), self._chunk_size
                )
            ]
        except:
            print("Supplied ID File has no ID's that can be extracted with the current settings.")
            print("Please loosen the settings and/or include more ID's in the supplied file.")
            exit(2)

        for chunk in self._id_chunks:
            keep_cols = [
                'Precursor m/z',
                'Identification Charge',
                'rt',
                'n_isos',
                # 'cf'
            ]
            chunk.drop(
                chunk.columns.difference(keep_cols),
                axis=1,
                inplace=True
            )
            chunk.reset_index(inplace=True)
            chunk.rename(
                columns={'index': 'id_index',
                         'Precursor m/z': 'mz',
                         'Identification Charge': 'z'},
                inplace=True
            )

    def _get_mzml_bounds(self, filename):
        '''Loads the mzml file (or at least a reader)

        Parameters
        ----------
        filename : str
            The name of the file containing the data

        Returns
        -------
        :obj:`bool`
            True if the file loads successfullly, False otherwise.

        Raises
        ------
        FileNotFoundError
            If the file does not exist

        '''
        # Open a pointer to the specified mzml file
        mzml_fp = pymzml.run.Reader(
            path_or_file=self.mzml_path,
            build_index_from_scratch=True
        )
        # Create a mapping of the relative index of each scan in the mzml
        #   to the native id given in the file. This allows us to iterate
        #   smoothly through the mzml file no matter how spaced out or
        #   oddly definied/formatted the given native indices are
        self._index_ID_map = {spec.index: spec.ID for spec in mzml_fp}
        mzml_bounds = dml.get_bounds(mzml_fp, self._index_ID_map)
        self._mzml_rt_min = mzml_bounds['rt_min']
        self._mzml_rt_max = mzml_bounds['rt_max']
        # make sure to close every file pointer that is not opened within a
        #   context manager (sometimes referred to as a 'with' block)
        mzml_fp.close()

    def run(self):
        '''Performs data extraction data from the mzml according to id file

        TODO: add a much longer explanation of why we needed special functions
              to chunk the data with overlap in the windows
        '''
        # These partial function applications (from python's 'functools'
        #   library) allow us to pass the same information to each of the
        #   chunks.
        func = partial(due.extract, self.settings_path)  # pass the settings
        func = partial(func, str(self.mzml_path))  # pass the mzml
        func = partial(func, str(self.id_path))  # pass the id path
        func = partial(func, self._index_ID_map)  # pass the index mapping

        if settings.debug_level == 0:
            try:
                results = list(
                    # tqdm creates our progress bar
                    tqdm(
                        self._mp_pool.imap_unordered(
                            func,
                            self._id_chunks
                        ),
                        total=len(self._id_chunks),
                        desc="Extracting: ",
                        leave=False
                    )
                )
            except Exception as e:
                print(
                    "The output for the data was too large to successfully output. Please use a smaller ID File to fix this issue.")

        if settings.debug_level >= 1:
            print('Beginning single-processor extraction.')
            results = []
            for chunk in tqdm(self._id_chunks, total=len(self._id_chunks), desc="Extracting: ", leave=False):
                results.append(func(chunk))

        # By filtering out the non-dataframe elements of this list, we exclude
        #   any incomplete or invalid chunks without ultimately terminating the
        #   program
        results = [r for r in results
                   if isinstance(r, pd.core.frame.DataFrame)]
        if len(results) == 0:
            # TODO: Implement and raise appropriate warning
            print('No data was successfully extracted')
        else:
            if len(results) == 1:
                self.model = one(results)
            else:
                self.model = pd.concat(results)
            drop_cols = ['rt', 'n_isos']
            self.ids = self.ids.drop(columns=drop_cols)
            self.model = self.ids.merge(
                self.model,
                how='outer',
                left_index=True,
                right_on='id_index'
            )

        self._mp_pool.close()
        self._mp_pool.join()


def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
