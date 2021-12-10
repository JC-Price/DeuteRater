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
        self.aa_labeling_dict = aa_label_df.loc[settings.study_type,].to_dict()

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
            self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
        #$breaks windows/python interactions if too many cores are used.  very niche application but still relevant
        if self._n_processors > 60:
            self.n_processors = 60
        self._mp_pool = mp.Pool(self._n_processors)
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
            func = partial(TheoryPreparer._mp_prepare, self.settings_path, aa_labeling_dict=self.aa_labeling_dict)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._enrichment_df), desc="Theory Generation: "
                )
            )
            
            
        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            results = []
            for row in tqdm(self._enrichment_df.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.Filename, sep='\t')
                
                func = partial(TheoryPreparer._mp_prepare, self.settings_path,
                                   aa_labeling_dict=self.aa_labeling_dict)
                if "mzs_list" in df.columns:
                    df.drop(inplace=True, columns=["mzs_list", "intensities_list", "rt_list", "baseline_list"])
                df = TheoryPreparer.func(df)
                df['time'] = row.Time
                df['enrichment'] = row.Enrichment
                df["sample_group"] = row.Sample_Group
                df["bio_rep"] = row.Biological_Replicate
                results.append(df)

        self.model = pd.concat(results)
        #if self.biomolecule_type == "Peptide":
        #    self.model = self.model.drop(columns=['drop'])
            
        self._mp_pool.close()
        self._mp_pool.join()

    @staticmethod
    def _mp_prepare(settings_path, args, aa_labeling_dict=None):
        settings.load(settings_path)
        #file_path, time, enrichment = args
        file_path, time, enrichment, sample_group, biological_replicate = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        if "mzs_list" in df.columns:
            df.drop(inplace=True, columns=["mzs_list", "intensities_list", "rt_list", "baseline_list"])
        df = TheoryPreparer._apply_filters(df)
        if aa_labeling_dict:
            #$don't include an else for either if statement.  no need to calculate if column exists
            #$ and we don't want to add the column if we can't calculate it since checking for it is an error check for later steps
            if literature_n_name not in df.columns:
                if aa_labeling_dict != "":
                    df = df.apply(TheoryPreparer._calculate_literature_n, axis =1 , args = (aa_labeling_dict,))
        df['time'] = time
        df['enrichment'] = enrichment
        df["sample_group"]  = sample_group
        df["bio_rep"] = biological_replicate
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

        #$get rid of lines that are missing mzs and abundances
        df = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        
        #$ we are going to drop all sequences that are too close in m/z and rt to other sequences
        #$swap to list of lists for speed
        df.sort_values(by='mz', inplace=True)
        mz_index = list(df.columns).index("mz") 
        rt_index = list(df.columns).index("rt")
        seq_index = list(df.columns).index(sequence_column_name)
        list_of_lists = df.values.tolist()
        too_close = []
        for i in range(len(list_of_lists)):
            for j in range(i+1, len(list_of_lists)):
                current_ppm = TheoryPreparer._ppm_calculator(list_of_lists[i][mz_index], list_of_lists[j][mz_index])
                if current_ppm > settings.mz_proximity_tolerance: 
                    break
                if abs(list_of_lists[i][rt_index] - list_of_lists[j][rt_index]) < settings.rt_proximity_tolerance and \
                        list_of_lists[i][seq_index] != list_of_lists[j][seq_index]:
                    too_close.extend([i,j])
        too_close = list(set(too_close))
        return df.drop(df.index[too_close])
    
    #$basic ppm calculator
    @staticmethod
    def _ppm_calculator(target, actual):
        ppm = (target-actual)/target * 1000000
        return abs(ppm)

def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
