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
import utils.NValueCalculator as nvct

import deuteconvert.peptide_utils as peputils

literature_n_name = "literature_n"
sequence_column_name = "Sequence"


class TheoryPreparer:
    def __init__(self, enrichment_path, out_path, settings_path, biomolecule_type):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        self.biomolecule_type = biomolecule_type
        if self.biomolecule_type == "Peptide":
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
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.75)
        else:
            self._n_processors = settings.n_processors
        # $breaks windows/python interactions if too many cores are used.  very niche application but still relevant
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
        results = []
        if settings.debug_level == 0:
            args_list = self._enrichment_df.to_records(index=False).tolist()
            if self.biomolecule_type == "Peptide":
                func = partial(TheoryPreparer._mp_prepare, self.settings_path, aa_labeling_dict=self.aa_labeling_dict)
            else:
                func = partial(TheoryPreparer._mp_prepare, self.settings_path)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._enrichment_df), desc="Theory Generation: "
                )
            )

        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            args_list = self._enrichment_df.to_records(index=False).tolist()
            for row in tqdm(self._enrichment_df.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.Filename, sep='\t')

                if self.biomolecule_type == "Peptide":
                    func = partial(TheoryPreparer._mp_prepare, self.settings_path,
                                   aa_labeling_dict=self.aa_labeling_dict)
                else:
                    func = partial(TheoryPreparer._mp_prepare, self.settings_path)
                if "mzs_list" in df.columns:
                    df.drop(inplace=True, columns=["mzs_list", "intensities_list", "rt_list", "baseline_list"])

                df = func(args_list[row.Index])
                df['timepoint'] = row.Time
                df['enrichment'] = row.Enrichment
                df["sample_group"] = row.Sample_Group
                df["bio_rep"] = row.Biological_Replicate
                results.append(df)

        self.model = pd.concat(results)
        # if self.biomolecule_type == "Peptide":
        #    self.model = self.model.drop(columns=['drop'])

        if settings.use_empir_n_value:
            self.model = self.model.reset_index(drop=True)
            self.model["row_num"] = np.arange(0, self.model.shape[0])
            self.model = self.model.loc[self.model["no_fn"] == ""]

            column_list = list(
                self.model.columns[
                    self.model.columns.isin(["Adduct", "sample_group", "Lipid Unique Identifier", "Sequence"])])
            column_list.sort()
            self.model["adduct_molecule_sg"] = self.model[column_list].agg("_".join, axis=1)

            # n_val_df = self.model
            calculator = nvct.NValueCalculator(self.model, self.settings_path, self.biomolecule_type)
            calculator.run()
            self.model = calculator.full_df

            full_df = self.model.copy()

            # Determine what the highest timepoint is and only look at those rows
            highest_timepoint = max(full_df['time'].unique())
            lipid_groups = full_df.groupby(by='adduct_molecule_sg')

            # # Compare reproducibility across reps
            for group in lipid_groups:
                group_df = group[1]
                high_tp_df = group_df.loc[group_df[
                                              'time'] == highest_timepoint]  # Only look at lipids that occur at the highest timepoint overall in the dataset. ie. D16 if timepoints are 0, 1, 4, 16
                if high_tp_df.empty:
                    # $ BN -1 is only for max time had no n-values (or grouping had no max time)
                    full_df.loc[full_df['adduct_molecule_sg'] == group[
                        0], 'n_value'] = -1  # If there is no lipids in the highest timepoint, set n_value as -1
                    continue
                # Remove reproducibility filter - CQ 15 Sept 2021
                if settings.remove_filters:
                    full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'n_value'] = round(
                        high_tp_df['empir_n'].median())
                else:
                    median_n = round(high_tp_df['empir_n'].median())  # BN rounding
                    # CQ Changed arrange so it has integers in the range. Trying to include as many values as possible within a range.
                    try:
                        median_range = np.arange(int(median_n - median_n * .1), round(median_n + median_n * .1) + 1,
                                                 1.0)  # BN swapped to a range added ", 1.0"
                    except:
                        pass
                    is_in_range_n = high_tp_df['empir_n'].apply(lambda x: x in median_range)
                    if is_in_range_n.all() and high_tp_df.shape[0] > 1:
                        all_n_values = list(high_tp_df['empir_n'])
                        if len(all_n_values) == 2:
                            all_n_values.append(np.median(all_n_values))
                        import scipy.stats as s
                        m, se = np.mean(all_n_values), s.sem(all_n_values)
                        if se == 0.0:
                            confidence_interval = (m, m)
                        else:
                            confidence_interval = s.t.interval(alpha=.90, df=len(all_n_values) - 1, loc=m, scale=se)

                        full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'n_value'] = median_n
                        full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'low_CI_n_value'] = confidence_interval[
                            0]
                        full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'high_CI_n_value'] = confidence_interval[
                            1]
                    elif high_tp_df.shape[0] == 1:
                        # If there is not 2 replicates of a specific lipid in the highest time course, set n_value as -2
                        full_df.loc[full_df['adduct_molecule_sg'] == group[
                            0], 'n_value'] = -2  # $ BN -2 indicates an error where max time n-values fell outside the "good"range
                    else:
                        # If the replicates of a specific lipid do not have reproducible n-values, set n_value as -3
                        full_df.loc[full_df['adduct_molecule_sg'] == group[0], 'n_value'] = -3

            full_df = full_df.rename(columns={'empir_n': 'n_val_calc_n',
                                              'n_value': 'empir_n'})

            full_df.loc[full_df.index, "n_val_calc_n"] = full_df["n_val_calc_n"]
            full_df.loc[full_df.index, "empir_n"] = full_df["empir_n"]
            self.model = full_df

        self._mp_pool.close()
        self._mp_pool.join()

    @staticmethod
    def _mp_prepare(settings_path, args, aa_labeling_dict=None):
        settings.load(settings_path)
        # file_path, time, enrichment = args
        file_path, time, enrichment, sample_group, biological_replicate = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        if "mzs_list" in df.columns:
            df.drop(inplace=True, columns=["mzs_list", "intensities_list", "rt_list", "baseline_list"])
        df = TheoryPreparer._apply_filters(df)
        if aa_labeling_dict:
            # $don't include an else for either if statement.  no need to calculate if column exists
            # $ and we don't want to add the column if we can't calculate it since checking for it is an error check for later steps
            if literature_n_name not in df.columns:
                if aa_labeling_dict != "":
                    df = df.apply(TheoryPreparer._calculate_literature_n, axis=1, args=(aa_labeling_dict,))
        df['time'] = time
        df['enrichment'] = enrichment
        df["sample_group"] = sample_group
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
        # This is 'clean_up_data' in the old deuterater
        df["no_fn"] = ""

        df.loc[df["mzs"].isna(), "no_fn"] = "no spectra extracted"

        data = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        data['drop'] = "False"
        # remove any rows that have an m/z that is within the proximity tolerance
        # Remove proximity Filter - CQ 15 Sept 2021
        if not settings.remove_filters:
            for row in data.itertuples():
                mask = ((data['mz'] - row.mz).abs() <
                        settings.mz_proximity_tolerance)
                data.loc[mask, 'drop'] = "True"
            # data = data[~data['drop']]

            df.loc[data.loc[data["drop"] == "True"].index] = "mz_proximity_tolerance_exceeded"

        # TODO: Check to see if no data went through

        return df


def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
