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

from functools import partial, reduce
import pandas as pd
import numpy as np
from scipy.stats import t  # pearsonr?
from scipy.optimize import curve_fit  # least_sq?
import warnings as w
import multiprocessing as mp
import concurrent.futures as cf

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.rate_equations as dur
from utils.graphing_tools import graph_rate

'''
Calculation of turnover rates
'''

group_column = "sample_group"


class RateCalculator:

    def __init__(self, model_path, out_path, graph_folders, settings_path, biomolecule_type):
        settings.load(settings_path)
        if model_path[-4:] == ".tsv":
            # 'Theoretical Unlabeled Normalized Abundances': str,
            self.model = pd.read_csv(model_path, sep='\t', low_memory=False, dtype={'theory_unlabeled_abunds': str,
                                                                                    'theory_labeled_abunds': str,
                                                                                    'normalized_empirical_abundances': str,
                                                                                    'low_labeling_peaks': str,
                                                                                    'frac_new_abunds': str,
                                                                                    'frac_new_abunds_std_dev': float})
        elif model_path[-4:] == ".csv":
            self.model = pd.read_csv(model_path, low_memory=False, dtype={'theory_unlabeled_abunds': str,
                                                                          'theory_labeled_abunds': str,
                                                                          'normalized_empirical_abundances': str,
                                                                          'low_labeling_peaks': str,
                                                                          'frac_new_abunds': str,
                                                                          'frac_new_abunds_std_dev': float})
        else:  # should never trigger unless we are fiddling with the gui
            raise ValueError("invalid file extension")
        self.out_path = out_path
        self.rate_model = None
        self.datapoint_model = None  # CQ
        self.graph_folders = graph_folders
        self.settings_path = settings_path

        # get the number of cores we're using for multiprocessing
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.75)
            # self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
        # breaks windows/python interactions if too many cores are used. very niche application but still relevant
        if self._n_processors > 60:
            self._n_processors = 60
        # self._mp_pool = mp.Pool(self._n_processors)

        self.biomolecule_type = biomolecule_type  # CQ

    def write(self):
        self.rate_model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )
        self.datapoint_model.to_csv(  # CQ
            path_or_buf=self.out_path[:-4] + "_datapoints.tsv",
            sep='\t',
            index=False
        )

    # this function governs which types of calculations to do and which
    #  rate equation to use. Collects the results and puts them in self.rate_model
    def calculate(self):
        #  catch the optimize warnings and so on
        w.filterwarnings("error")
        rate_results = []
        datapoint_results = []  # CQ
        max_time = max(self.model["time"])

        # get the appropriate analyte column
        if self.biomolecule_type == "Peptide":
            # peptide_analyte_id_column is the Protein_ID column
            id_col = settings.peptide_analyte_id_column
        elif self.biomolecule_type == "Lipid":
            # lipid_analyte_id_column is the Lipid Unique Identifier column
            id_col = settings.lipid_analyte_id_column
        else:
            print("Unknown Analyte Type")
            return

        # could put manual bias correction here, would be slightly more
        # efficient, but make the logic less clear

        # Decide if rates should be put together with adducts or not.
        if self.biomolecule_type == "Lipid" and settings.separate_adducts:
            self.model[id_col] = self.model[id_col] + self.model["Adduct"].apply(lambda x: f"_{x}")

        # group_column is the sample_group column
        groups = self.model.groupby(by=[id_col, group_column])

        # Prepare for fit
        if settings.asymptote == 'Variable':
            p0 = [.1, 1]
            rate_eq = partial(
                dur.simple,
                p_adj=settings.proliferation_adjustment
            )
        elif settings.asymptote == 'Fixed':
            p0 = .1
            rate_eq = partial(
                dur.simple,
                p_adj=settings.proliferation_adjustment,
                a=settings.fixed_asymptote_value
            )
        else:
            p0 = np.nan
            rate_eq = np.nan

        # we're going to load up a function for mp.  we need to load up the
        # function. don't add fn_col, fn_std_dev, calc_type, and manual_bias,
        # they will be added later
        rate_function = partial(RateCalculator._mp_function,
                                settings_path=self.settings_path,
                                graph_folders=self.graph_folders,
                                rate_eq=rate_eq,
                                max_time=max_time,
                                p0=p0)

        # call the rate for each relevant measurement type
        # add measurement specific arguments to partial function

        # perform calculations for an abundance based rate
        # settings.debug_level = 1
        if settings.fraction_new_calculation == "abundance" or settings.fraction_new_calculation == "combined":
            temp_rate_function = partial(rate_function,
                                         calc_type="Abundance",
                                         fn_col="abund_fn",
                                         fn_std_dev="frac_new_abunds_std_dev",
                                         manual_bias=settings.abundance_manual_bias,
                                         std_dev_filter=settings.abundance_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                    results = list(tqdm(executor.map(temp_rate_function, groups), total=len(groups),
                                        desc="Abundance Rate Calculation: ", leave=True))
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Abundance Rate Calculation: "):
                    results.append(temp_rate_function(group))

            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ

        # perform calculations for a neutromer spacing based rate
        if settings.fraction_new_calculation == "neutromer spacing" or settings.fraction_new_calculation == "combined":
            temp_rate_function = partial(rate_function,
                                         calc_type="Spacing",
                                         fn_col="nsfn",
                                         fn_std_dev="frac_new_mzs_std_dev",
                                         manual_bias=settings.spacing_manual_bias,
                                         std_dev_filter=settings.spacing_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                    results = list(tqdm(executor.map(temp_rate_function, groups), total=len(groups),
                                        desc="Spacing Rate Calculation: ", leave=True))
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Spacing Rate Calculation: "):
                    results.append(temp_rate_function(group))
            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ

        # perform calculations for a combined rate (neutromer spacing + abundance)
        if settings.fraction_new_calculation == "combined":
            temp_rate_function = partial(rate_function,
                                         calc_type="Combined",
                                         fn_col="cfn",
                                         fn_std_dev="frac_new_combined_std_dev",
                                         manual_bias=settings.combined_manual_bias,
                                         std_dev_filter=settings.combined_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                    results = list(tqdm(executor.map(temp_rate_function, groups), total=len(groups),
                                        desc="Combined Rate Calculation: ", leave=True))
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Combined Rate Calculation: "):
                    results.append(temp_rate_function(group))
            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ

        # from https://stackoverflow.com/questions/44327999/python
        # -pandas-merge-multiple-dataframes answer 1 accessed 9/16/2020
        # merges any number of dfs into one df
        self.rate_model = reduce(lambda left, right: pd.merge(left, right,
                                                              on=['analyte_id', 'analyte_name', "group_name"],
                                                              how='outer'), rate_results)
        self.datapoint_model = pd.concat(datapoint_results[0])
        # swap back to normal warning behaviour
        w.filterwarnings("default")

        if not settings.verbose_rate:
            self.trim_verbose_data()

    def trim_verbose_data(self):
        needed_columns = ['analyte_id', 'analyte_name', "group_name"]
        potential_columns = ["{} rate", "{} 95pct_confidence", "{} half life"]
        if settings.asymptote != "Fixed":
            potential_columns.insert(1, "{} asymptote")
        if settings.fraction_new_calculation == "abundance" or settings.fraction_new_calculation == "combined":
            for p in potential_columns:
                needed_columns.append(p.format("Abundance"))
        if settings.fraction_new_calculation == "neutromer spacing" or settings.fraction_new_calculation == "combined":
            for p in potential_columns:
                needed_columns.append(p.format("Spacing"))
        if settings.fraction_new_calculation == "combined":
            for p in potential_columns:
                needed_columns.append(p.format("Combined"))
        self.rate_model = self.rate_model[needed_columns]

    # function passed into the multiprocessing for calculations
    def _mp_function(data_tuple, settings_path, fn_col, fn_std_dev, calc_type, manual_bias, std_dev_filter,
                     graph_folders,
                     rate_eq, max_time, p0, biomolecule_type):
        w.filterwarnings("error")
        settings.load(settings_path)
        pd.options.mode.chained_assignment = None

        # id_values will have the Protein_ID and sample_group
        # group will be the dataframe containing the actual data
        id_values, group = data_tuple[0], data_tuple[1]
        id_name = id_values[0]
        sample_group_name = id_values[1]

        # grab the common name for the group (based on analyte name column)
        if biomolecule_type == "Peptide":
            common_name = group[settings.peptide_analyte_name_column].iloc[0]
        elif biomolecule_type == "Lipid":
            common_name = group[settings.lipid_analyte_name_column].iloc[0]
        else:
            raise Exception("Not a Valid Biomolecule Type")

        # drop error string could do earlier for more speed, but this is
        # clearer and allows errors that affect only one calculation type
        group = RateCalculator._error_trimmer(group, [fn_col, fn_std_dev])

        # the copy is just to avoid a SettingWithCopy warning in a few
        # operations.  if it causes problems remove and suppress warning
        if not settings.remove_filters:
            group = group[group[fn_std_dev] < std_dev_filter].copy()

        # return error message if there is no data for group
        if len(group) == 0:
            result = RateCalculator._make_error_message(
                "No Isotope Envelopes Agree", "", id_name, common_name,
                sample_group_name, calc_type, 0, 0, 0, 0)
            return result, group

        if settings.use_empir_n_value:
            # check if cv is above the set filter value
            if group.iloc[0]['cv'] == "no valid time points" or group.iloc[0]['cv'] == "error occurred":
                result = RateCalculator._make_error_message("Error: see n_value column", "",
                                                            id_name, common_name, sample_group_name, calc_type,
                                                            0, 0, 0, 0)
                return result, group
            elif float(group.iloc[0]['cv']) > settings.n_value_cv_limit:
                result = RateCalculator._make_error_message(f"cv is greater than {settings.n_value_cv_limit}", "",
                                                            id_name, common_name, sample_group_name, calc_type,
                                                            0, 0, 0, 0)
                return result, group

        # make sure we have the right dtypes for important columns
        try:
            convert_dict = {'n_value': float, 'time': float}
            if settings.fraction_new_calculation == "abundance" or settings.fraction_new_calculation == "combined":
                convert_dict['abund_fn'] = float
                convert_dict['frac_new_abunds_std_dev'] = float
            if settings.fraction_new_calculation == "neutromer spacing" or settings.fraction_new_calculation == "combined":
                convert_dict['nsfn'] = float
                convert_dict['frac_new_mzs_std_dev'] = float
            if settings.fraction_new_calculation == "combined":
                convert_dict['cfn'] = float
                convert_dict['frac_new_combined_std_dev'] = float
            group = group.astype(convert_dict)
        except Exception as e:
            result = RateCalculator._make_error_message(
                "value could not be determined", "Fraction new value(s) could not be used for rate calculations",
                id_name, common_name, sample_group_name, calc_type, np.NaN,
                np.NaN, np.NaN, np.NaN)
            return result, group

        # offset all values by a certain amount (instrument bias)
        if settings.bias_calculation == "calculated":
            bias = RateCalculator._calc_bias(group, fn_col)
            group[fn_col] = group[fn_col] - bias
        elif settings.bias_calculation == "manual":  # user designated bias
            group[fn_col] = group[fn_col] - manual_bias

        # set time as the x-axis and fraction new as the y-axis
        xs = np.concatenate(([0], group['time'].to_numpy()))
        ys = np.concatenate(([settings.y_intercept_of_fit], group[fn_col].to_numpy()))

        # if settings.use_outlier_removal:
        #     # ys = RateCalculator.mask_outliers(ys)
        #     xs, ys, devs = RateCalculator._roll(xs, ys)
        # else:
        #     devs = np.concatenate(([settings.error_of_zero], group[fn_std_dev].to_numpy()))

        # remove any outliers if we have more than 4 timepoints
        if settings.use_outlier_removal and biomolecule_type == "Lipid":
            ys, indices = RateCalculator.mask_outliers(ys)
            xs = xs[indices]
        elif biomolecule_type == "Peptide" and settings.use_outlier_removal:
            # xs, ys, devs = RateCalculator._roll(xs, ys)
            ys, indices = RateCalculator.mask_outliers(ys)
            xs = xs[indices]

        # Get the number of unique time points, and continue if not enough
        num_unique_times = len(set(group['time']))

        if biomolecule_type == "Peptide":
            unique_length = len(set(group[settings.unique_sequence_column]))
            # for lipids or metabolites or similar, the unique length means
            # nothing.  if needed can add the elif
        else:
            unique_length = ""
        num_measurements = len(group.index)

        # TODO: this is not ideal but this is a good first attempt
        num_files = len(set(group["mzml_path"]))

        # TODO: Handle Num Bio Reps in Stuff

        # I think this fixes the issues with technical replicates.
        # need to use astype(str) on all or if someone uses numbers for group names or replicate names issues result
        num_bio_reps = len(
            set(group["time"].astype(str) + group["sample_group"].astype(str) + group["bio_rep"].astype(str)))

        if num_unique_times < settings.minimum_nonzero_points:
            result = RateCalculator._make_error_message(
                "Insufficient times", "", id_name, common_name, sample_group_name,
                calc_type, num_measurements, num_unique_times, unique_length, num_files)
            return result, group

        # perform fit
        try:
            # DO NOT use std dev as the Sigma because it creates influential outliers
            # don't use sigma unless we have a different
            # popt are the optimal values for the parameters so the sum of the squared residuals of f(xdata, *popt) - ydata is minimized
            # pcov is the estimated approximate covariance of popt.
            # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html for more details
            try:
                # Having some issues with curve_fit. Should we try Lmfit or the lm option?
                # https://stackoverflow.com/questions/50371428/scipy-curve-fit-raises-optimizewarning-covariance-of-the-parameters-could-not
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
                popt, pcov = curve_fit(f=rate_eq, xdata=xs, ydata=ys, p0=p0)
            except Exception as e:
                if type(e) == ValueError:
                    current_exception = "Incompatible x and/or y data"
                elif type(e) == RuntimeError:
                    current_exception = "least_squares minimization failed"
                else:
                    current_exception = "covariance of parameters couldn't be estimated"

                result = RateCalculator._make_error_message(
                    "value could not be determined", current_exception, id_name,
                    common_name, sample_group_name, calc_type, num_measurements,
                    num_unique_times, unique_length, num_files)

                return result, group

            # pull results of fit into variables
            rate = popt[0]
            asymptote = popt[1] if len(popt) > 1 else settings.fixed_asymptote_value
            # TODO ci uses degrees of freedom = n-k where n is the number of points and k is the number of parameters estimated
            # including intercept in linear regression.  if asymptote is fixed k=1 otherwise k=2 (intercept is fit by equation, not data)
            # not counting charge states and different peptides as unique measurements.
            # despite the claim in the documentation, according to statistics consultation and every site I checked, the np.sqrt(np.diag(pcov))[0]
            # is standard error, not std dev, so don't divide by sqrt of n

            # try to find confidence interval. Sometimes covariance will be negative, in which case we won't be able to calculate one.
            # If so, we'll set confint to np.NaN and we won't draw any error lines
            confint = None
            for value in np.diag(pcov):
                if value < 0:
                    confint = np.NaN
                    break
            if confint is None:
                try:
                    confint = t.ppf(.975, num_files - len(popt)) * np.sqrt(np.diag(pcov))[0]
                except Exception as e:
                    confint = np.NaN

            try:
                y_predicted = dur.simple(xs, rate, asymptote, settings.proliferation_adjustment)
                r_2 = dur.calculate_r2(ys, y_predicted)
            except Exception as e:
                result = RateCalculator._make_error_message(
                    "r_2 calculation failed", e, id_name,
                    common_name, sample_group_name, calc_type, num_measurements,
                    num_unique_times, unique_length, num_files)
                return result, group

            result = {
                'analyte_id': id_name,
                'analyte_name': common_name,
                'group_name': sample_group_name,
                '{} rate'.format(calc_type): rate,
                '{} asymptote'.format(calc_type): asymptote,
                '{} 95pct_confidence'.format(calc_type): confint,
                '{} half life'.format(calc_type): RateCalculator._halflife(rate),
                '{} R2'.format(calc_type): r_2,
                "{} files observed in".format(calc_type): num_files,
                '{} num_measurements'.format(calc_type): num_measurements,
                '{} num_time_points'.format(calc_type): num_unique_times,
                '{} uniques'.format(calc_type): unique_length,
                '{} exceptions'.format(calc_type): "",
                "rate_graph_time_points_x": "; ".join(map(str, list(xs))),
                "normed_isotope_data_y": "; ".join(map(str, list(ys)))
            }

            if confint is not np.NaN:
                result['{} std_error'.format(calc_type)] = np.sqrt(np.diag(pcov))[0]
            else:
                result['{} std_error'.format(calc_type)] = np.NaN

            # if there is an asymptote need to provide it
            if biomolecule_type == "Peptide":
                graph_name = "{}_{}_{}".format(id_name, sample_group_name,
                                               fn_col)
                graph_title = "{}_{}_{}\nk={}, a={}, r-squared={}, ".format(common_name,
                                                                            sample_group_name, fn_col,
                                                                            result[f'{calc_type} rate'],
                                                                            1.0, round(r_2, 3))
            elif biomolecule_type == "Lipid":
                graph_name = "{}_{}_{}".format(id_name[:-1],
                                               sample_group_name, fn_col)
                graph_title = "{}_{}_{}\nk={:.3f}, a={:.3f}, r-squared={}".format(id_name[:-1],
                                                                                  sample_group_name, fn_col,
                                                                                  result[f'{calc_type} rate'],
                                                                                  result['{} rate'.format(calc_type)],
                                                                                  round(r_2, 3))
                if "dietary_lit_n_avg" in list(group.columns):
                    lit_n_str = ""
                    if group.iloc[0]["dietary_lit_n_avg"] != 0:
                        diet_lit_n_range = [float(a) for a in group.iloc[0]["dietary_lit_n_range"][1:-1].split(", ")]
                        lit_n_str += "({:.1f}-{:.1f}) (dietary), ".format(diet_lit_n_range[0], diet_lit_n_range[1])
                    if group.iloc[0]["de_novo_lit_n_avg"] != 0:
                        de_novo_lit_n_range = [float(a) for a in group.iloc[0]["de_novo_lit_n_range"][1:-1].split(", ")]
                        lit_n_str += "({:.1f}-{:.1f}) (de_novo)".format(de_novo_lit_n_range[0], de_novo_lit_n_range[1])
                    if lit_n_str != "":
                        graph_title += "\nlit-n={}".format(lit_n_str)
            else:
                graph_name = None
                graph_title = None
            try:
                if r_2 > 0.5:
                    graph_folder = graph_folders[0]
                else:
                    graph_folder = graph_folders[1]
                if settings.graph_output_format == "none":
                    print("No graphs were generated. Change DeuteRater settings if you want graphs to be generated.")
                elif biomolecule_type == "Lipid":
                    graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                               rate_eq, graph_folder, max_time,
                               settings.asymptote, calc_type=fn_col, full_data=group, title=graph_title)
                else:
                    # if biomolecule_type == "Peptide":
                    #     graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                    #                rate_eq, graph_folder, max_time,
                    #                settings.asymptote, devs, title=graph_title)
                    # else:
                    graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                               rate_eq, graph_folder, max_time,
                               settings.asymptote, title=graph_title)
            except Exception as c:
                print("Graphing failed with the following error: ")
                print(c)
        except Exception as c:
            # "we have a guess but are unsure" warning
            if type(c).__name__ == "OptimizeWarning":
                current_exception = 'OptimizeWarning: optimal fit could not be found'
            # couldn't find the minimum
            elif type(c).__name__ == "RuntimeError":
                current_exception = 'fit could not be found'
            else:
                raise c  # will stop here so don't need to consider further
            result = RateCalculator._make_error_message(
                "value could not be determined", current_exception, id_name,
                common_name, sample_group_name, calc_type, num_measurements,
                num_unique_times, unique_length, num_files)

        return result, group

    # function that does a mad outlier check on a numpy array
    @staticmethod
    def mask_outliers(y_values):
        if len(y_values) == 1:
            return y_values, [0]
        med = np.median(y_values)
        differences = abs(y_values - med)
        med_abs_dev = np.median(differences)
        # print("outlier data: ")
        # print(y_values)
        # print(med)
        # print(differences)
        # print(med_abs_dev)
        z_values = (differences * .6745) / med_abs_dev
        # good_indices = np.where(z_values < settings.zscore_cutoff)[0]
        good_indices = np.asarray(z_values < settings.zscore_cutoff).nonzero()[0]
        return y_values[good_indices], good_indices

    # rollup time points
    @staticmethod
    def _roll(old_x, old_y):
        new_x = np.unique(old_x)
        new_y, errors = [], []
        for x in new_x:
            #  grab indices for the time
            # x_time_indices = np.where(old_x == x)[0]
            x_time_indices = np.asarray(old_x == x).nonzero()[0]

            # grab only relevant indices of y and outlier check
            if len(old_y) >= 4:
                temp_y_values = RateCalculator.mask_outliers(old_y[x_time_indices])

                # second length check (first is in _mask_outliers) to avoid
                # medians and std devs on statistically invalid data
                if len(temp_y_values) == 1:
                    new_y.extend(temp_y_values)
                    errors.append(settings.error_of_non_replicated_point)
                else:
                    new_y.append(np.median(temp_y_values))
                    standard_dev = np.std(temp_y_values, ddof=1)
                    # I have seen std dev be 0 which causes problems for the fit, so we'll just force no 0s
                    if standard_dev == 0.0:
                        errors.append(settings.error_of_zero)
                    else:
                        errors.append(standard_dev)
                return new_x, np.array(new_y), np.array(errors)
            else:
                return old_x, old_y, []

    # need to calculate median value at time zero in the target column
    # already checked for std dev error at this point
    @staticmethod
    def _calc_bias(df, target_column):
        test_df = df[df["time"] == 0]
        if len(test_df) > 0:
            return np.median(test_df[target_column])
        else:
            return 0

    # implemented and working
    @staticmethod
    def _halflife(x):
        return np.log(2) / x

    # just to shorten the code for the error messages
    @staticmethod
    def _make_error_message(main_error, exception, id_name, common_name, group_name,
                            calc_type, num_m, num_times, unique_length, num_files):
        result = {'analyte_id': id_name,
                  'analyte_name': common_name,
                  'group_name': group_name,
                  '{} rate'.format(calc_type): main_error,
                  '{} asymptote'.format(calc_type): main_error,
                  '{} std_error'.format(calc_type): main_error,
                  '{} 95pct_confidence'.format(calc_type): main_error,
                  '{} half life'.format(calc_type): main_error,
                  '{} R2'.format(calc_type): main_error,
                  "{} files observed in".format(calc_type): num_files,
                  '{} num_measurements'.format(calc_type): num_m,
                  '{} num_time_points'.format(calc_type): num_times,
                  '{} uniques'.format(calc_type): unique_length,
                  '{} exceptions'.format(calc_type): exception
                  }
        return result

    # removes errors from a previous analysis step, which are strings
    # just cut out those rows df is the dataframe to trim, error columns is a list of columns where
    # errors are located. returns the trimmed dataframe
    @staticmethod
    def _error_trimmer(df, error_columns):
        df[error_columns] = df[error_columns].apply(pd.to_numeric, errors='coerce')
        temp = df.dropna(axis=0, subset=error_columns)
        return temp
