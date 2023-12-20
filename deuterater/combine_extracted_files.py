# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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
After the data is extracted form mzmls we need to prepare it for future analysis
At this point the user has provided information on which files go with which
experimental subject and the enrichment, as well as which extracted filtes to consider

Therefore this will add the users data on subjects and deuterium enrichment,
do some basic filtering, calculate a n value based on aa_labeling_sites.tsv
(in the resources folder) and then merge all into one file

may merge with initial_intensity_calculator, which has a similar basic filtering
role and occurs immediately afterwards
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


#$as with all the calculation steps this is a class for consistent calls in the main
class CombineExtractedFiles():
    def __init__(self, enrichment_path, out_path, settings_path, needed_columns):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        self.needed_columns = needed_columns
        #$collect necessary components to determine n_value from amino acids
        aa_label_df = pd.read_csv(settings.aa_labeling_sites_path, sep='\t')
        aa_label_df.set_index('study_type', inplace=True)
        self.aa_labeling_dict = aa_label_df.loc[settings.label_key, ].to_dict()
        
        
        #$pull in two sub tables from the output table. posibilities for .tsv and .csv files
        if self.enrichment_path.suffix == '.tsv':
            self._file_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t',
                usecols = ["Filename", "Time", "Subject ID"]
            )
            self._enrichment_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t',
                usecols = ["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
            )
        elif self.enrichment_path.suffix == '.csv':
            self._file_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=',',
                usecols = ["Filename", "Time", "Subject ID"]
            )
            self._enrichment_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=',',
                usecols = ["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
            )
            
        #$since the multiple sub tables can have different length, get rid
        #$of the rows that are empty
        self._file_data.dropna(inplace = True, how = "all")
        self._enrichment_data.dropna(inplace = True, how = "all")
        self._data_dict = self.collect_enrichment_data()
        
        #$if multiprocessing need to set that up. more than 60 cores causes problems for windows
        if settings.recognize_available_cores is True:
            self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
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
    #$read in data from the user
    def collect_enrichment_data(self):
        data_dict = {}
        for subject, subject_df in self._enrichment_data.groupby("Subject ID Enrichment"):
            x_values = ", ".join([str(x) for x in subject_df["Time Enrichment"]])
            y_values = ", ".join([str(y) for y in subject_df["Enrichment"]])
            data_dict[subject] = [x_values, y_values]
        return data_dict
            
    #$run the analysis.  this function doesn't have any calculation itself (other than merging results)
    #$it prepares a function for multiprocessing and thne begins the multiprocessing
    def prepare(self):
        if settings.debug_level == 0:
            args_list = self._file_data.to_records(index=False).tolist()
            func = partial(CombineExtractedFiles._mp_prepare, self.settings_path, self._data_dict, self.aa_labeling_dict)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._file_data),
                    desc="Combine Extracted Files: "
                )
            
            )
            
            
        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            results = []
            for row in tqdm(self._file_data.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.file, sep='\t')
                df = CombineExtractedFiles._apply_filters(df)
                if literature_n_name not in df.columns:
                    if self.aa_labeling_dict != "":
                        df = df.apply(CombineExtractedFiles._calculate_literature_n, axis =1 , args = (self.aa_labeling_dict,))
                
                df['time'] = row.time
                df["sample_id"]  = row.sample_id
                df["Time Enrichment"] = self._data_dict[row.sample_id][0]
                df["Enrichment Values"] = self._data_dict[row.sample_id][1]
                results.append(df)

        self.model = pd.concat(results)
        #$now we need to filter columns
        #$otherwise the carry forward increases file size quite a bit
        #$by doing here it should not affect anything.
        self.model = self.model[self.needed_columns]
        

        self._mp_pool.close()
        self._mp_pool.join()

    #$actually runs the relevant calculation. Yes reloading the settings is necessary
    #$because each process has its own global variables in windows
    @staticmethod
    def _mp_prepare(settings_path, data_dict,  aa_labeling_dict, args):
        settings.load(settings_path)
        #file_path, time, enrichment = args
        file_path, time, sample_id = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        df = CombineExtractedFiles._apply_filters(df)
        #$ if the user or a previous process defined n, that's fine.  but it will be
        #$needed in the next step so calculate if necessary.
        if literature_n_name not in df.columns:
            if aa_labeling_dict != "":
                df = df.apply(CombineExtractedFiles._calculate_literature_n, axis =1 , args = (aa_labeling_dict,))
        df['time'] = time
        df["sample_id"]  = sample_id
        df["Time Enrichment"] = data_dict[sample_id][0]
        df["Enrichment Values"] = data_dict[sample_id][1]
            
        return df
    
    #$calculate the n value based on amino acid sequence
    @staticmethod
    def _calculate_literature_n(row, aa_labeling_dict):
        aa_counts = {}
        for aa in row["Sequence"]:
            if aa not in aa_counts.keys():
                aa_counts[aa] = 0
            aa_counts[aa] += 1
        literature_n = peputils.calc_add_n(aa_counts, aa_labeling_dict)
        row[literature_n_name] = literature_n
        return row

    #$because the extractor may assign different ids to the same peak in different files,
    #$we need to account for that. if two peaks are too clos in both retention time and mz
    #$remove both as we can't be sure of the id without ms/ms which we aren't dealing with
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
        #$take out of the dataframe for speed
        df.sort_values(by='mz', inplace=True)
        mz_index = list(df.columns).index("mz") 
        rt_index = list(df.columns).index("rt")
        seq_index = list(df.columns).index("Sequence")
        list_of_lists = df.values.tolist()
        too_close = []
        for i in range(len(list_of_lists)):
            for j in range(i+1, len(list_of_lists)):
                current_ppm = CombineExtractedFiles._ppm_calculator(list_of_lists[i][mz_index], list_of_lists[j][mz_index])
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
    
#$can't really use as a main 
def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
