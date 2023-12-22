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

"""
as with combine_extracted_files, this provides some basic filtering and 
calculations that are relatively simple on their own, but need to be done and it's
easier here than in the main cacluation function

The main filters here are n value and sequence length must be sufficiently large
to be confident we are going to actually see good signal to noise. Finally we must
use the emass algorightm to calculate the unenriched relative isotopic envelope

may merge into combine_extracted_files.py
"""



from tqdm import tqdm
import pandas as pd
import multiprocessing as mp
import numpy as np

from pathlib import Path
from functools import partial

from utils.emass import emass
import deuterater.settings as settings


max_isos = 5 # constant based on the n_isos based on the mass (done in the extractor)
p0_guess = 1 # seems to work for most fits. if it causes problems we can adjust


# as with all the calculation steps this is a class for consistent calls in the main
class theoretical_enrichment_calculator(object):
    def __init__(self, prepared_data_path, out_path, settings_path):
        settings.load(settings_path)
        self.settings_path = settings_path
        
        self.prepared_data_path = Path(prepared_data_path)
        self.out_path = out_path
        
        if self.prepared_data_path.suffix == '.tsv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep='\t'
            )
        elif self.prepared_data_path.suffix == '.csv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep=','
            )
            # if multiprocessing need to set that up. more than 60 cores causes problems for windows
            if settings.recognize_available_cores is True:
                # BD: Issue with mp.cpu_count() finding too many cores available
                self._n_processors = round(mp.cpu_count() * 0.80)
            else:
                self._n_processors = settings.n_processors
            if self._n_processors > 60:
                self.n_processors = 60
            
    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )    
        
    # run the analysis.  this function doesn't have any calculation itself (other than merging, transposing, and setting index of the results)
    # it prepares a function for multiprocessing and thne begins the multiprocessing
    def prepare(self):
        unique_sequnces_df = self.data_df.drop_duplicates(subset = ["Sequence"])
    
        new_columns = theoretical_enrichment_calculator._make_new_columns()
        func = partial(theoretical_enrichment_calculator._individual_process, 
                       new_columns = new_columns,
                       minimum_n_value = settings.min_allowed_n_values,
                       minimum_sequence_length = settings.min_aa_sequence_length)
        df_split = np.array_split(unique_sequnces_df, len(unique_sequnces_df))
                     
        mp_pools = mp.Pool(self._n_processors)
        
        final_df = pd.concat(tqdm(
            mp_pools.imap(func, df_split), 
            total = len(df_split),
            desc="Theory Generation: "
            ),axis =1)
        
        mp_pools.close()
        mp_pools.join()
        final_df = final_df.T
        final_df = final_df.set_index("Sequence")
        self.model = pd.merge(self.data_df, final_df, left_on= "Sequence", right_index = True)
        
    # actually runs the relevant calculation. 
    @staticmethod
    def _individual_process(df, new_columns, 
                            minimum_n_value,minimum_sequence_length):
         variable_list = []
         for row in df.itertuples():
            output_series = pd.Series(index = new_columns, dtype = "object")
            output_series["Sequence"] = row.Sequence
            # drop rows we will not use before doing any comples calculations
            if len(row.Sequence) < minimum_sequence_length:
                variable_list.append(theoretical_enrichment_calculator._error_message_results(
                    f"Sequence is less than {minimum_sequence_length} amino acids",
                    output_series))
                continue
            if row.literature_n < minimum_n_value:
                variable_list.append(theoretical_enrichment_calculator._error_message_results(
                    f"less than {minimum_n_value} labeling sites",
                    output_series))
                continue
            intensity_values = \
                theoretical_enrichment_calculator._fit_emass(row.cf,
                      row.n_isos
                )

            output_series["Theoretical Unlabeled Normalized Abundances"] = ", ".join(intensity_values)
            variable_list.append(output_series)
         return(pd.concat(variable_list,axis =1))
     

    # if an error happens it is most efficient to have a function
    def _error_message_results(error_message, output_series):
        # don't need to know which names are which or how many columns there are, 
        # just need python to fill all non-Sequence columns
        # position 0 is the sequence name which we don't wish to overwrite
        for index_name in output_series.index[1:]:
            output_series[index_name] = error_message
        return output_series

                
    

    # calculate unlabeled intensity, if we need to return  m/z values or
    # adjust for different n_values, do it here or in emass itself.
    def _fit_emass(sequence, n_isos):
        intensity_values = emass(
                    sequence,
                    n_isos
                )
        return [str(i) for i in intensity_values]

    # this creates the header for the variables. is a function in case we need 
    # to add various columns  (if we want to graph emass output or something)
    @staticmethod
    def _make_new_columns():
        new_columns = ["Sequence"]
        new_columns.extend(["Theoretical Unlabeled Normalized Abundances"])
        return new_columns
    