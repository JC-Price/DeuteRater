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
import utils.NValueCalculator as nvct


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
        #TODO: is this for peptides only or do we use this same process for lipids?
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

        # Brought this over from another DR version. Used to calculate n-values when we aren't using literature ones.
        # It uses the NValueCalculator.py and n_value_calc_emass.py files in the Util folder - Ben Driggs
        #TODO: should this be based on molecule type or use_empir_n_value setting?
        if settings.use_empir_n_value:
            self.model = self.model.reset_index(drop=True)
            self.model["row_num"] = np.arange(0, self.model.shape[0])
            self.model = self.model.loc[self.model["no_fn"] == ""]

            column_list = list(
                self.model.columns[
                    self.model.columns.isin(["Adduct", "sample_group", "Lipid Unique Identifier", "Sequence"])])
            column_list.sort()
            self.model["adduct_molecule_sg"] = self.model[column_list].agg("_".join, axis=1)

            # We don't want to calculate N-values for Day 0 data or for enrichment less than 0.005 - Ben Driggs
            # n_val_df = self.model
            calculator = nvct.NValueCalculator(self.model, self.settings_path, self.biomolecule_type)
            calculator.run()
            self.model = calculator.full_df

            full_df = self.model.copy()

            # Determine what the highest time point is and only look at those rows
            highest_timepoint = max(full_df['time'].unique())
            lipid_groups = full_df.groupby(by='adduct_molecule_sg')

            # TODO: make sure we are calculating the median from all timepoints the user put "yes" in calculate_n_value column

            # # Compare reproducibility across reps
            for group in lipid_groups:
                group_df = group[1]
                high_tp_df = group_df.loc[group_df[
                                              'time'] == highest_timepoint]  # Only look at lipids that occur at the highest time point overall in the dataset. ie. D16 if timepoints are 0, 1, 4, 16

                # Find median n-value from all timepoints where the user put "yes" in calculate_n_value column
                calc_n_value_df = group_df.loc[group_df['calculate_n_value'] == "yes"]

                if calc_n_value_df.empty:
                    # $ BN -1 is only for max time had no n-values (or grouping had no max time)
                    full_df.loc[full_df['adduct_molecule_sg'] == group[
                        0], 'n_value'] = -1  # If there is no lipids in the highest timepoint, set n_value as -1
                    continue
                # Remove reproducibility filter - CQ 15 Sept 2021
                if settings.remove_filters:
                    full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'n_value'] = round(
                        calc_n_value_df['empir_n'].median())
                else:
                    median_n = round(calc_n_value_df['empir_n'].median())  # BN rounding
                    # CQ Changed arrange so that it has integers in the range. Trying to include as many values as possible within a range.
                    try:
                        median_range = np.arange(int(median_n - median_n * .1), round(median_n + median_n * .1) + 1,
                                                 1.0)  # BN swapped to a range added ", 1.0"
                    except:
                        pass
                    is_in_range_n = calc_n_value_df['empir_n'].apply(lambda x: x in median_range)
                    if is_in_range_n.all() and calc_n_value_df.shape[0] > 1:
                        all_n_values = list(calc_n_value_df['empir_n'])
                        if len(all_n_values) == 2:
                            all_n_values.append(np.median(all_n_values))
                        import scipy.stats as s
                        m, se = np.mean(all_n_values), s.sem(all_n_values)
                        if se == 0.0:
                            confidence_interval = (m, m)
                        else:
                            confidence_interval = s.t.interval(alpha=.90, df=len(all_n_values) - 1, loc=m, scale=se)

                        full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'n_value'] = median_n
                        # .loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = median_n
                        # full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'low_CI_n_value'] = confidence_interval[
                        #     0]
                        # full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'high_CI_n_value'] = confidence_interval[
                        #     1]
                        full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'low_CI_n_value'] = confidence_interval[0]
                        full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'high_CI_n_value'] = confidence_interval[1]
                    elif calc_n_value_df.shape[0] == 1:
                        # If there is not 2 replicates of a specific lipid in the highest time course, set n_value as -2
                        full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = -2  # $ BN -2 indicates an error where max time n-values fell outside the "good"range
                    else:
                        # If the replicates of a specific lipid do not have reproducible n-values, set n_value as -3
                        full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = -3

            full_df = full_df.rename(columns={'empir_n': 'n_val_calc_n',
                                              'n_value': 'empir_n'})

            full_df.loc[full_df.index, "n_val_calc_n"] = full_df["n_val_calc_n"]
            full_df.loc[full_df.index, "empir_n"] = full_df["empir_n"]
            self.model = full_df
        
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
    