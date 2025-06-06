# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, Benjamin Driggs,
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
import numpy as np  # noqa
from copy import copy
import multiprocessing as mp
import concurrent.futures as cf

from functools import partial

from tqdm import tqdm  # noqa: 401

from utils.emass import parse_cf, fn_emass, normalize

import deuterater.settings as settings


# at this point we have far too many columns to be practical
# this is not a problem for theory since that is largely a combination of
# the table and the extracted files, so it is a nice place to store that data
# at this point we just need the bare essentials for calculation to ease
# troubleshooting and visibility
columns_to_drop = ["Precursor Retention Time (sec)", "rt_start", "rt_end",
                   "rt_width", "Precursor m/z", "Identification Charge",
                   "id_index", "lookback_mzs", "lookback_abundances",
                   "lookahead_mzs", "lookahead_abundances", "rt_min",
                   "rt_max", "mads",
                   ]

peptide_extra_columns_to_drop = ["quality", "avg_ppm", "start_loc",
                                 "end_loc", "species", 'gene_name',
                                 "protein_existence", "sequence_version"]

# itertuples doesn't like some characters
# these characters include () and space
# this is because the elements are called using periods
# therefore we will declare any particularly problematic.
# we'll force that to change here
# if we're swapping between literature_n and theory_n use this dictionary
protein_itertuple_renamer = {
    "Protein ID": "Protein_ID",
    "Homologous Proteins": "Homologous_Proteins",
    "Protein Name": "Protein_Name",
}


class FractionNewCalculator:
    def __init__(self, model_path, out_path, settings_path, biomolecule_type):
        settings.load(settings_path)
        self.settings_path = settings_path
        if biomolecule_type == "Peptide":
            itertuple_renamer = copy(protein_itertuple_renamer)
            # settings.use_empir_n_value = False
        else:
            itertuple_renamer = None

        self.biomolecule_type = biomolecule_type

        self.error = ""

        # get the number of cores we're using for multiprocessing
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.75)
            # self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
        # breaks windows/python interactions if too many cores are used. V-ery niche application but still relevant
        if self._n_processors > 60:
            self._n_processors = 60

        # self._mp_pool = mp.Pool(self._n_processors)
        
        if model_path[-4:] == ".tsv":
            self.model = pd.read_csv(model_path, sep='\t', low_memory=False)
        elif model_path[-4:] == ".csv":
            self.model = pd.read_csv(model_path, low_memory=False)
        else:  # should never trigger unless we are fiddling with the gui
            raise ValueError("invalid file extension")

        self.out_path = out_path
        if itertuple_renamer:
            self.model = self.model.rename(columns=itertuple_renamer)

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False)

    def generate(self):
        # just drop the unneeded columns
        self.model = self.model.drop(columns_to_drop, axis=1, errors="ignore")

        if "no_fn" in self.model.columns:
            self.model = self.model.loc[pd.isna(self.model["no_fn"])]
        if self.biomolecule_type == "Peptide":
            self.model = self.model.drop(peptide_extra_columns_to_drop, axis=1,
                                         errors="ignore")
        if self.model["n_value"].isnull().all().all():
            self.error = ("N value is not present.  "
                          "Please ensure the n value has been provided in the "
                          "literature_n column or has been calculated in the "
                          "n_value column. If the proper column is filled "
                          "please check your settings for using the proper "
                          "column")
            return

        if settings.fraction_new_calculation == "abundance" or settings.fraction_new_calculation == "combined":
            self.model = self.model.assign(
                theory_unlabeled_abunds="",
                theory_labeled_abunds="",
                normalized_empirical_abundances="",
                low_labeling_peaks="",
                frac_new_abunds="",
                frac_new_abunds_std_dev="",
                abund_fn=""
            )
        if settings.fraction_new_calculation == "neutromer spacing" or settings.fraction_new_calculation == "combined":
            self.model = self.model.assign(
                observed_neutral_masses="",
                theory_unlabeled_mzs="",
                theory_labeled_mzs="",
                frac_new_mzs="",
                frac_new_mzs_outlier_checked="",
                frac_new_mzs_std_dev="",
                nsfn=""
            )
        if settings.fraction_new_calculation == "combined":
            self.model = self.model.assign(
                frac_new_combined="",
                frac_new_combined_outlier_checked="",
                frac_new_combined_std_dev="",
                cfn=""
            )

        # now that we have dropped useless columns we can start the multiprocessing
        # don't use np.split, or it will error if chunks are not equal sized
        model_pieces = np.array_split(self.model, len(self.model))

        if settings.debug_level == 0:
            func = partial(FractionNewCalculator._mp_prepare)
            func = partial(func, settings_path=self.settings_path)
            func = partial(func, biomolecule_type=self.biomolecule_type)

            with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                results = list(
                    tqdm(executor.map(func, model_pieces), total=len(self.model),
                         desc="Fraction New Calculation: ", leave=True))

        if settings.debug_level >= 1:
            print('Beginning single-processor fraction new calculation.')
            results = []
            for row in tqdm(model_pieces, total=len(self.model), desc="Fraction New Calculation: "):
                results.append(FractionNewCalculator._mp_prepare(row, self.settings_path, self.biomolecule_type))

        self.model = pd.concat(results, sort=False)

    # generic error output for filters based on n values or other criteria
    # preferably will do before the calculations so can leave all blank except
    #  actual measurement
    @staticmethod
    def _error_method(df, row, error_message):
        if settings.fraction_new_calculation == "abundance":
            df.at[row.Index, "abund_fn"] = error_message
        if settings.fraction_new_calculation == "neutromer spacing":
            df.at[row.Index, "nsfn"] = error_message
        if settings.fraction_new_calculation == "combined":
            df.at[row.Index, "cfn"] = error_message
        return df

    @staticmethod
    def _mp_prepare(df, settings_path, biomolecule_type):
        settings.load(settings_path)

        # Remove NaN values for row.n_value and replace with error code -4
        df['n_value'] = df['n_value'].fillna(int(-4))

        # can start with itertuples. If need be can swap to apply
        for row in df.itertuples(index=True):

            if settings.use_empir_n_value:
                # first check if the cv is above the allowed amount
                if row.cv == "no valid time points" or row.cv == "error occurred":
                    df = FractionNewCalculator._error_method(df, row, "Error: see n_value column")
                    continue

            if settings.remove_filters:
                if row.n_value < 0:
                    df = FractionNewCalculator._error_method(df, row, "N value is less than {}".format(settings.min_allowed_n_values))
                    continue
            else:
                # Do any initial filtering to save time calculating
                if float(row.n_value) < settings.min_allowed_n_values:
                    df = FractionNewCalculator._error_method(df, row, "N value is less than {}".format(settings.min_allowed_n_values))
                    continue
            if biomolecule_type == "Peptide" and len(row.Sequence) < settings.min_aa_sequence_length:
                df = FractionNewCalculator._error_method(df, row, "Fewer than {} amino acids".format(settings.min_aa_sequence_length))
                continue

            # if the user chooses 0 as enrichment it will cause divide by zero
            # error later, so we will use the settings to force it
            # result should be close to 0 change anyway
            if float(row.enrichment) != 0.0:
                use_enrich = float(row.enrichment)
            else:
                # Usually is 0.05 so no need to store it as a setting
                # use_enrich = settings.enrichment_of_zero
                use_enrich = 0.05

            try:
                num_h, parsed_cf = parse_cf(row.cf_w_adduct)
            except:
                num_h, parsed_cf = parse_cf(row.cf)

            # generate emass information for the row
            _, enriched_results = fn_emass(
                parsed_cf=parsed_cf,
                n_list=[0, float(row.n_value)],
                n_H=num_h,
                low_pct=0,
                high_pct=use_enrich,
                num_peaks=int(row.n_isos),
                testing=False
            )
            e_mzs, e_abunds = enriched_results

            if settings.fraction_new_calculation == "abundance" or settings.fraction_new_calculation == "combined":
                # takes the e_abunds values and puts them into a single string value
                FractionNewCalculator._prepare_row(df, row, "abunds", e_abunds)

                # normalizes the peaks by ratio using the normalize function in emass.py
                normalized_empirical_abunds = FractionNewCalculator._normalize_abundances(row.abundances[1:-1].split(", "))

                # stores normalized abundances to dataframe
                df.at[row.Index, 'normalized_empirical_abundances'] = normalized_empirical_abunds

                # calculate deltas for theory abundances
                theory_abund_deltas = FractionNewCalculator._calculate_deltas(e_abunds.loc[1][1:], e_abunds.loc[0][1:])

                # calculate deltas for empirical abundances
                empirical_abund_deltas = FractionNewCalculator._calculate_deltas(
                    [float(x) for x in normalized_empirical_abunds.split(", ")],
                    e_abunds.loc[0][1:]
                )

                # takes the theory and empirical abundance deltas and removes any peaks that are below 0.04
                theory_abund_deltas, empirical_abund_deltas, removed_peaks = \
                    FractionNewCalculator._trim_abunds(theory_abund_deltas, empirical_abund_deltas)

                # store what peaks were removed in the low_labeling_peaks list
                df.at[row.Index, "low_labeling_peaks"] = removed_peaks

                # TODO: Should this be included?
                # df.at[row.Index, 'delta_I'] = np.std(empirical_abund_deltas)

                # don't need to break if we only have one or zero. combined and
                # spacing are still fine, we just shouldn't do the other calculations here
                if len(theory_abund_deltas) < 2:
                    minimum_abund_change = 0.04
                    df.at[row.Index, "abund_fn"] = f"Insufficient peaks with theory above {minimum_abund_change}"
                    all_frac_new_abunds = []
                else:
                    # calculate fraction based off empirical and theory abundance deltas
                    all_frac_new_abunds = FractionNewCalculator._calculate_fractions(
                        empirical_abund_deltas, theory_abund_deltas
                    )

                    df = FractionNewCalculator.final_calculations(df, row, "abunds", all_frac_new_abunds,
                                                                  "abund_fn", theory_abund_deltas[0])

            if settings.fraction_new_calculation == "neutromer spacing" or settings.fraction_new_calculation == "combined":
                FractionNewCalculator._prepare_row(df, row, "mzs", e_mzs)
                theory_mz_deltas = FractionNewCalculator._calculate_deltas(
                    FractionNewCalculator._mz_deltas(e_mzs.loc[1][1:]),
                    FractionNewCalculator._mz_deltas(e_mzs.loc[0][1:])
                )
                # need to trim the parenthesis from the row.mzs string
                # also turn it into floats to do math on it
                observed_mzs = [float(x) for x in row.mzs[1:-1].split(", ")]
                observed_neutral_mass = [y * int(row.z)
                                         for y in observed_mzs]
                df.at[row.Index, "observed_neutral_masses"] = ", ".join([str(x) for x in observed_neutral_mass])
                empirical_mz_deltas = FractionNewCalculator._calculate_deltas(
                    FractionNewCalculator._mz_deltas(observed_neutral_mass),
                    FractionNewCalculator._mz_deltas(e_mzs.loc[0][1:])
                )
                all_frac_new_mzs = FractionNewCalculator._calculate_fractions(empirical_mz_deltas, theory_mz_deltas)
                df = FractionNewCalculator.final_calculations(df, row, "mzs", all_frac_new_mzs, "nsfn")

            if settings.fraction_new_calculation == "combined":
                all_frac_new_combined = all_frac_new_abunds + all_frac_new_mzs
                df = FractionNewCalculator.final_calculations(df, row, "combined", all_frac_new_combined, "cfn")

            # filter out any rows that have a fraction new standard deviation greater than the max_fn_standard deviation setting
            if settings.fraction_new_calculation == "abundance":
                df = FractionNewCalculator.stddev_filter(df, 'frac_new_abunds_std_dev', 'abund_fn')
            elif settings.fraction_new_calculation == "neutromer spacing":
                df = FractionNewCalculator.stddev_filter(df, 'frac_new_mzs_std_dev', 'nsfn')
            else:
                df = FractionNewCalculator.stddev_filter(df, "frac_new_combined_std_dev", 'cfn')

        return df

    @staticmethod
    def stddev_filter(df, fn_stddev_column, fn_column):
        # for any rows that have an n_value_stddev greater than the max_fn_standard_deviation, remove fraction new and stddev
        if isinstance(df.iloc[0][fn_column], float):
            if float(df.iloc[0][fn_stddev_column]) > settings.max_fn_standard_deviation:
                if fn_column == "abund_fn":
                    # df.frac_new_abunds = ""
                    # df.frac_new_abunds_std_dev = ""
                    df.abund_fn = f"{fn_stddev_column} greater than {settings.max_fn_standard_deviation}"
                if fn_column == "nsfn":
                    # df.frac_new_mzs = ""
                    # df.frac_new_mzs_std_dev = ""
                    df.nsfn = f"{fn_stddev_column} greater than {settings.max_fn_standard_deviation}"
                if fn_column == "cfn":
                    # df.frac_new_abunds = ""
                    # df.frac_new_abunds_std_dev = ""
                    # df.abund_fn = ""
                    # df.frac_new_mzs = ""
                    # df.frac_new_mzs_std_dev = ""
                    # df.nsfn = ""
                    # df.frac_new_combined = ""
                    # df.frac_new_combined_std_dev = ""
                    df.cfn = f"{fn_stddev_column} greater than {settings.max_fn_standard_deviation}"
        return df

    @staticmethod
    def _series_to_string(pd_series):
        string_values = [str(a) for a in pd_series.values]
        return " ".join(string_values)

    @staticmethod
    def _mz_deltas(A):
        return [mz_value - A[0] for mz_value in A[1:]]

    @staticmethod
    def _calculate_deltas(A, B):
        return [a - b for a, b in zip(A, B)]

    # calculate fraction based off empirical and theory abundance deltas
    # if fraction new > 1, it means empirical delta is > than 1
    # if fraction new < 0, it means empirical delta is negative
    @staticmethod
    def _calculate_fractions(empirical_deltas, theory_deltas):
        # observed over expected
        return [a / b for a, b in zip(empirical_deltas, theory_deltas)]

    # normalizes the peaks by ratio using the normalize function in emass.py
    @staticmethod
    def _normalize_abundances(abundance_list):
        normalized_list = normalize([float(a) for a in abundance_list])
        return ", ".join([str(n) for n in normalized_list])

    # need to remove the abundances with a low maximum. usually caused by the
    # peak rising and falling. since we are calculating based on unidirectional
    # change this can cause problems for the std dev filter
    @staticmethod
    def _trim_abunds(theory, empirical):
        drop_list, new_theory, new_empirical = [], [], []
        # TODO ask JC if this is something we would need changed very often - Ben D
        minimum_abund_change = 0.04
        for t in range(len(theory)):
            if abs(theory[t]) < minimum_abund_change:
                drop_list.append(f'M{t}')
            else:
                new_theory.append(theory[t])
                new_empirical.append(empirical[t])
        return new_theory, new_empirical, ", ".join(drop_list)

    # for right now this is not necessary, but we may need to expand it later
    #  can add a length check or outlier check to this later
    @staticmethod
    def _std_dev_calculator(list_to_std_dev):
        return np.std(list_to_std_dev, ddof=1)

    # need to add an outlier check to the neutromer spacing (mzs)
    # as well as any combined data.  This is a median absolute
    @staticmethod
    def _mad_outlier_check(data_to_check, z_score_cutoff):
        initial_median = np.median(data_to_check)  # only calculate once
        differences = [abs(x - initial_median) for x in data_to_check]
        mad = np.median(differences)
        zscores = [.6745 * d / mad for d in differences]
        # even if we could do a list comprehension for the next step
        # it would be too long for easy readability
        good_turnovers = []
        for z in range(len(zscores)):
            if zscores[z] < z_score_cutoff:
                good_turnovers.append(data_to_check[z])
        if len(good_turnovers) > 1:
            return good_turnovers
        else:
            return data_to_check

    # initial columns used for abundance and spacing.  currently equations are
    # too different to put in but may do later
    @staticmethod
    def _prepare_row(full_df, row, keyword, df):
        # need to do [1:] in order to cut out the n_D
        full_df.at[row.Index, 'theory_unlabeled_{}'.format(keyword)] = FractionNewCalculator._series_to_string(df.loc[0][1:])
        full_df.at[row.Index, 'theory_labeled_{}'.format(keyword)] = FractionNewCalculator._series_to_string(df.loc[1][1:])
        return full_df

    # compresses the code to make the last few output columns
    @staticmethod
    def final_calculations(df, row, keyword, data, final_column, m0_max_delta=1):
        df.at[row.Index, 'frac_new_{}'.format(keyword)] = ", ".join([str(r) for r in data])

        # NOTE: data = all_frac_new_abunds
        if keyword == "abunds":
            min_allowed_abund_max_delta = 0.04
            df.at[row.Index, 'frac_new_{}_std_dev'.format(keyword)] = FractionNewCalculator._std_dev_calculator(data)
            if abs(m0_max_delta) < min_allowed_abund_max_delta:
                df.at[row.Index, final_column] = "Low Theoretical M0 delta possible. Point rejected"
            else:
                if settings.abundance_type == "M0":
                    df.at[row.Index, final_column] = data[0]
                elif settings.abundance_type == "Average":
                    df.at[row.Index, final_column] = np.median(data)
                elif settings.abundance_type == "Highest":
                    try:
                        abundances = [float(a) for a in row.abundances[1:-1].split(", ")]
                        abunds_to_drop = df.low_labeling_peaks.iloc[0]
                        # handle a single drop
                        if len(abunds_to_drop) == 2:
                            abunds_to_drop = int(abunds_to_drop[1])
                            abundances.pop(abunds_to_drop)
                        elif len(abunds_to_drop) == 0:
                            pass
                        else:
                            # TODO: what does this mean? - Ben D
                            print("Need to implement multiple drops... help!")

                        df.at[row.Index, final_column] = data[np.argmax(abundances)]
                    except:
                        print("Need to implement multiple drops... help!")
        else:
            # TODO: need to do further analysis on outlier check for these
            data = FractionNewCalculator._mad_outlier_check(data, settings.zscore_cutoff)
            df.at[row.Index, 'frac_new_{}_outlier_checked'.format(keyword)] = ", ".join([str(r) for r in data])
            df.at[row.Index, 'frac_new_{}_std_dev'.format(keyword)] = FractionNewCalculator._std_dev_calculator(data)
            df.at[row.Index, final_column] = np.median(data)
        return df
