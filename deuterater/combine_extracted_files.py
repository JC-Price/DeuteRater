# -*- coding: utf-8 -*-
"""
Copyright (c) 2024 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, Benjamin Driggs,
    Coleman Nielsen, J.C. Price, and Brigham Young University
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

import traceback  # noqa: 401
from pathlib import Path
from functools import partial
import multiprocessing as mp
import concurrent.futures as cf

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.NValueCalculator as nvct
import utils.experimental_n_value_calculator as envct

import resources.peptide_utils as peputils

'''
After the data is extracted form mzmls we need to prepare it for future analysis
At this point the user has provided information on which files go with which
experimental subject and the enrichment, as well as which extracted filters to consider

Therefore this will add the users data on subjects and deuterium enrichment,
do some basic filtering, calculate a literature n value based on aa_labeling_sites.tsv (for peptides) or 
an empirical n value based on theory using NValueCalculator.py (for lipids) and then merge them all into one file

may merge with initial_intensity_calculator, which has a similar basic filtering
role and occurs immediately afterwards
'''

literature_n_name = "literature_n"


# as with all the calculation steps this is a class for consistent calls in the main
class CombineExtractedFiles:
    def __init__(self, enrichment_path, out_path, settings_path, needed_columns, biomolecule_type, graph_folder):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        self.needed_columns = needed_columns
        self.biomolecule_type = biomolecule_type
        self.graph_folder = graph_folder

        # settings.use_empir_n_value = True

        # collect necessary components to determine n_value from amino acids
        if self.biomolecule_type == "Peptide":
            aa_label_df = pd.read_csv(settings.aa_labeling_sites_path, sep='\t')
            aa_label_df.set_index('study_type', inplace=True)
            self.aa_labeling_dict = aa_label_df.loc[settings.label_key,].to_dict()
            # settings.use_empir_n_value = False
        
        # since the multiple sub tables can have different length, get rid
        # of the rows that are empty
        # Get data from time and enrichment table
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

        # self._file_data.dropna(inplace=True, how="all")
        # self._enrichment_data.dropna(inplace=True, how="all")
        # self._data_dict = self.collect_enrichment_data()
        
        # if multiprocessing need to set that up. more than 60 cores causes problems for windows
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.80)
        else:
            self._n_processors = settings.n_processors
        if self._n_processors > 60:
            self._n_processors = 60

        # self._mp_pool = mp.Pool(self._n_processors)
        self.out_path = out_path
        self.model = None

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    # read in data from the user
    def collect_enrichment_data(self):
        data_dict = {}
        for subject, subject_df in self._enrichment_data.groupby("Subject ID Enrichment"):
            x_values = ", ".join([str(x) for x in subject_df["Time Enrichment"]])
            y_values = ", ".join([str(y) for y in subject_df["Enrichment"]])
            data_dict[subject] = [x_values, y_values]
        return data_dict
        
    # run the analysis.  this function doesn't have any calculation itself (other than merging results)
    # it prepares a function for multiprocessing and then begins the multiprocessing

    def prepare(self):
        results = []
        if settings.debug_level == 0:
            self._enrichment_df['biomolecule_type'] = self.biomolecule_type
            args_list = self._enrichment_df.to_records(index=False).tolist()
            if self.biomolecule_type == "Peptide":
                func = partial(CombineExtractedFiles._mp_prepare, self.settings_path, aa_labeling_dict=self.aa_labeling_dict)
            else:
                func = partial(CombineExtractedFiles._mp_prepare, self.settings_path)

            with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                results = list(
                    tqdm(executor.map(func, args_list), total=len(args_list), desc="Theory Generation: ",
                         leave=True))

        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            self._enrichment_df['biomolecule_type'] = self.biomolecule_type
            args_list = self._enrichment_df.to_records(index=False).tolist()
            for row in tqdm(self._enrichment_df.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.Filename, sep='\t')

                if self.biomolecule_type == "Peptide":
                    func = partial(CombineExtractedFiles._mp_prepare, self.settings_path,
                                   aa_labeling_dict=self.aa_labeling_dict)
                else:
                    func = partial(CombineExtractedFiles._mp_prepare, self.settings_path)
                if "mzs_list" in df.columns:
                    df.drop(inplace=True, columns=["mzs_list", "intensities_list", "rt_list", "baseline_list"])

                df = func(args_list[row.Index])
                df['timepoint'] = row.Time
                df['enrichment'] = row.Enrichment
                df["sample_group"] = row.Sample_Group
                df["bio_rep"] = row.Biological_Replicate
                df["calculate_n_value"] = row.Calculate_N_Value
                results.append(df)

        self.model = pd.concat(results)
        
        # if there is a column already named n_value, we want to remove it. Otherwise, when we merge the columns from the n-value
        # calculations to the model, it will append _x to the original n_value column name and _y to the n_value_calculator n_value column.
        if self.model.columns.isin(['n_value']).any():
            self.model.drop(columns=['n_value'], inplace=True)

        if settings.use_empir_n_value:
            self.model = self.model.reset_index(drop=True)
            self.model["row_num"] = np.arange(0, self.model.shape[0])

            # Used this for the human version to exclude rows that had really low intensities
            # self.model = self.model.loc[self.model["no_fn"] == ""]
            if self.biomolecule_type == 'Peptide':
                column_list = list(self.model.columns[self.model.columns.isin(
                    ["sample_group", "Sequence", 'Protein ID'])])
                column_list.sort()
                self.model["molecule_sg"] = self.model[column_list].agg("_".join, axis=1)
            else:
                column_list = list(self.model.columns[self.model.columns.isin(["Adduct", "sample_group", "Lipid Unique Identifier", "Sequence", 'Protein ID'])])
                column_list.sort()
                self.model["adduct_molecule_sg"] = self.model[column_list].agg("_".join, axis=1)

            # We don't want to calculate N-values for Day 0 data or for enrichment less than 0.005 - Ben Driggs
            # n_val_df = self.model
            calculator = nvct.NValueCalculator(self.model, self.settings_path, self.biomolecule_type, self.out_path, graphs_location=self.graph_folder)
            # calculator = envct.Experimental_NValueCalculator(self.model, self.settings_path, self.biomolecule_type, self.out_path, graphs_location=self.graph_folder)
            calculator.run()
            self.model = calculator.full_df

            full_df = self.model.copy()

            # Determine what the highest time point is and only look at those rows
            # highest_timepoint = max(full_df['time'].unique())
            # lipid_groups = full_df.groupby(by='adduct_molecule_sg')

            # TODO: Since we already have calculated the average n-value, what aspects of this reproducibility code do we need? - Ben D
            # # Compare reproducibility across reps
            # for group in lipid_groups:
            #     group_df = group[1]
            #
            #     # TODO: Do we remove reproducibility filter - CQ 15 Sept 2021
            #     if settings.remove_filters:
            #         full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'n_value'] = round(
            #             group_df['empir_n'].median())
            #     else:
            #         median_n = round(group_df['empir_n'].median())  # BN rounding
            #         # CQ Changed arrange so that it has integers in the range. Trying to include as many values as possible within a range.
            #         try:
            #             median_range = np.arange(int(median_n - median_n * .1), round(median_n + median_n * .1) + 1,
            #                                      1.0)  # BN swapped to a range added ", 1.0"
            #         except:
            #             pass
            #         is_in_range_n = group_df['empir_n'].apply(lambda x: x in median_range)
            #         if is_in_range_n.all() and group_df.shape[0] > 1:
            #             all_n_values = list(group_df['empir_n'])
            #             if len(all_n_values) == 2:
            #                 all_n_values.append(np.median(all_n_values))
            #             import scipy.stats as s
            #             m, se = np.mean(all_n_values), s.sem(all_n_values)
            #             if se == 0.0:
            #                 confidence_interval = (m, m)
            #             else:
            #                 confidence_interval = s.t.interval(alpha=.90, df=len(all_n_values) - 1, loc=m, scale=se)
            #
            #             full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'n_value'] = median_n
            #             # .loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = median_n
            #             # full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'low_CI_n_value'] = confidence_interval[
            #             #     0]
            #             # full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'high_CI_n_value'] = confidence_interval[
            #             #     1]
            #             full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'low_CI_n_value'] = confidence_interval[0]
            #             full_df.loc[(full_df['adduct_molecule_sg'] == group[0]), 'high_CI_n_value'] = confidence_interval[1]
            #         elif group_df.shape[0] == 1:
            #             # If there is not 2 replicates of a specific lipid in the highest time course, set n_value as -2
            #             full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = -2  #  BN -2 indicates an error where max time n-values fell outside the "good"range
            #         else:
            #             # If the replicates of a specific lipid do not have reproducible n-values, set n_value as -3
            #             full_df.loc[(full_df['adduct_molecule_sg'] == group[0]) & (full_df['calculate_n_value'] == "yes"), 'n_value'] = -3

            # full_df = full_df.rename(columns={'empir_n': 'empir_n_value'})

            # full_df = full_df.rename(columns={'empir_n': 'n_value', 'n_value': 'empir_n'})

            # full_df.drop('empir_n', axis=1, inplace=True)

            # full_df.loc[full_df.index, "n_value"] = full_df["n_value"]

            full_df.drop("calculate_n_value", axis=1, inplace=True)
            self.model = full_df
        else:
            # put literature_n values into n_value column
            self.model['n_value'] = self.model['literature_n']

    # actually runs the relevant calculation. Yes reloading the settings is necessary
    # because each process has its own global variables in windows
    @staticmethod
    def _mp_prepare(settings_path, args, aa_labeling_dict=None):
        settings.load(settings_path)
        # file_path, time, enrichment = args
        # file_path, time, sample_id = args
        file_path, time, enrichment, sample_group, biological_replicate, calculate_n_values, biomolecule_type = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        df = CombineExtractedFiles._apply_filters(df)
        #  if the user or a previous process defined n, that's fine.  but it will be
        # needed in the next step so calculate if necessary.
        if biomolecule_type == "Peptide" and literature_n_name not in df.columns:
            if aa_labeling_dict != "":
                df = df.apply(CombineExtractedFiles._calculate_literature_n, axis=1, args=(aa_labeling_dict,))

        df['time'] = time
        df['enrichment'] = enrichment
        df["sample_group"] = sample_group
        df["bio_rep"] = biological_replicate
        df['calculate_n_value'] = calculate_n_values
        return df
    
    # calculate the n value based on amino acid sequence
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

    # because the extractor may assign different ids to the same peak in different files,
    # we need to account for that. if two peaks are too clos in both retention time and mz
    # remove both as we can't be sure of the id without ms/ms which we aren't dealing with
    @staticmethod
    def _apply_filters(df):
        """filters the internal dataframe

        This function does not modify the dataframe in place.

        Parameters
        ----------
        df : :obj:`pandas.Dataframe`
            The internal dataframe of the theory_value_prep function

        Returns
        -------
        :obj:`pandas.Dataframe`
            The filtered dataframe. Does not modify in place.
        """
        # get rid of lines that are missing mzs and abundances
        df = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        
        #  we are going to drop all sequences that are too close in m/z and rt to other sequences
        # take out of the dataframe for speed
        # df.sort_values(by='mz', inplace=True)
        # mz_index = list(df.columns).index("mz")
        # rt_index = list(df.columns).index("rt")
        # seq_index = list(df.columns).index("Sequence")
        # list_of_lists = df.values.tolist()
        # too_close = []
        # for i in range(len(list_of_lists)):
        #     for j in range(i+1, len(list_of_lists)):
        #         current_ppm = CombineExtractedFiles._ppm_calculator(list_of_lists[i][mz_index], list_of_lists[j][mz_index])
        #         if current_ppm > settings.mz_proximity_tolerance:
        #             break
        #         if abs(list_of_lists[i][rt_index] - list_of_lists[j][rt_index]) < settings.rt_proximity_tolerance and \
        #                 list_of_lists[i][seq_index] != list_of_lists[j][seq_index]:
        #             too_close.extend([i, j])
        # too_close = list(set(too_close))
        # return df.drop(df.index[too_close])
        return df
    
    # basic ppm calculator
    @staticmethod
    def _ppm_calculator(target, actual):
        ppm = (target-actual)/target * 1000000
        return abs(ppm)


# can't really use as a main 
def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
