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

'''Theoretical Value Preparation

This purpose of this module is to calculate the expected theoretical values
of each identified analyte.

'''

import pandas as pd
import numpy as np  # noqa: 401

import traceback  # noqa: 401
from pathlib import Path
from functools import partial
import multiprocessing as mp

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings

import deuteconvert.peptide_utils as peputils


literature_n_name = "literature_n"
sequence_column_name = "Sequence"

class TheoryPreparer():
    def __init__(self, enrichment_path, out_path, settings_path):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        
        aa_label_df = pd.read_csv(settings.aa_label_path, sep='\t')
        aa_label_df.set_index('study_type', inplace=True)
        self.aa_labeling_dict = aa_label_df.loc[settings.study_type, ].to_dict()
        
        if self.enrichment_path.suffix == '.tsv':
            self._enrichment_df = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t'
            )
        elif self.enrichment_path.suffix == '.csv':
            self._enrichment_df = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=','
            )
        if settings.recognize_available_cores is True:
            self._n_partitions = mp.cpu_count()
        else:
            self._n_partitions = settings.n_partitions

        self._mp_pool = mp.Pool(self._n_partitions)
        self.out_path = out_path
        self.model = None

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    def prepare(self):
        if settings.debug_level == 0:
            args_list = self._enrichment_df.to_records(index=False).tolist()
            func = partial(TheoryPreparer._mp_prepare, self.settings_path, self.aa_labeling_dict)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._enrichment_df)
                )
            
            )
            
            
        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            results = []
            for row in tqdm(self._enrichment_df.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.file, sep='\t')
                df = TheoryPreparer._apply_filters(df)
                #$don't include an else for either if statement.  no need to calculate if column exists
                #$ and we don't want to add the column if we can't calculate it since checking for it is an error check for later steps
                if literature_n_name not in df.columns:
                    df = df.apply(TheoryPreparer._calculate_literature_n, axis =1 , args = (self.aa_labeling_dict,))
                df['time'] = row.time
                df['enrichment'] = row.enrichment
                df["sample_group"]  = row.sample_group
                results.append(df)

        self.model = pd.concat(results)
        self.model = self.model.drop(columns=['drop'])

        self._mp_pool.close()
        self._mp_pool.join()

    @staticmethod
    def _mp_prepare(settings_path, aa_labeling_dict, args):
        settings.load(settings_path)
        #file_path, time, enrichment = args
        file_path, time, enrichment, sample_group = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        df = TheoryPreparer._apply_filters(df)
        #$don't include an else for either if statement.  no need to calculate if column exists
        #$ and we don't want to add the column if we can't calculate it since checking for it is an error check for later steps
        if literature_n_name not in df.columns:
            if aa_labeling_dict != "":
                df = df.apply(TheoryPreparer._calculate_literature_n, axis =1 , args = (aa_labeling_dict,))
        df['time'] = time
        df['enrichment'] = enrichment
        df["sample_group"]  = sample_group
        return df

    @staticmethod
    def _calculate_literature_n(row, aa_labeling_dict):
        aa_counts = {}
        for aa in row[sequence_column_name]:
            if aa not in aa_counts.keys():
                aa_counts[aa] = 0
            aa_counts[aa] += 1
        literature_n = peputils.calc_add_n(aa_counts, aa_labeling_dict)
        row[literature_n_name] = literature_n
        return row

    # def _load(self):
    #     '''Pulls in the relevant data from the model.

    #     Parameters
    #     ----------
    #     model : :obj:`pandas.Dataframe`
    #         The unmodified data model

    #     Returns
    #     -------
    #     :obj:`pandas.Dataframe`
    #         A dataframe containing a copy of the relevant columns

    #     '''
    #     # TODO: Basic checks like whether the data looks right need done
    #     if not isinstance(self.model, pd.DataFrame):
    #         try:
    #             # TODO: csv/tsv flexibility?
    #             self.model = pd.read_csv(self.model, sep='\t')
    #         except Exception as e:
    #             # TODO: better exception logging
    #             print(e)
    #             traceback.print_tb(e.__traceback__)
    #             raise
    #     try:
    #         # TODO: csv/tsv flexibility?
    #         self._enrichment_df = pd.read_csv(self.enrichment_path, sep='\t')
    #     except Exception as e:
    #         # TODO: better exception logging
    #         print(e)
    #         traceback.print_tb(e.__traceback__)
    #         raise

    @staticmethod
    def _apply_filters(df):
        '''filters the internal dataframe

        This function does not modify the dataframe in place.

        Parameters
        ----------
        df : :obj:`pandas.Dataframe`
            The internal dataframe of the theory_value_prep function

        Returns
        -------
        :obj:`pandas.Dataframe`
            The filtered dataframe. Does not modify in place.
        '''
        # This is 'clean_up_data' in the old deuterater
        # This is a
        data = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        data['drop'] = False
        for row in data.itertuples():
            mask = ((data['mz'] - row.mz).abs() <
                    settings.mz_proximity_tolerance)
            data.loc[mask, 'drop'] = True
        data = data[~data['drop']]

        # TODO: Check to see if no data went through
        return data


def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
