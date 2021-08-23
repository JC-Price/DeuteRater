# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor, J.C. Price, and Brigham Young University
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


'''Calculation of turnover rates

'''

from functools import partial, reduce
import pandas as pd
import numpy as np
from scipy.stats import t  # pearsonr?
from scipy.optimize import curve_fit  # least_sq?
import warnings as w
import multiprocessing as mp

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.rate_equations as dur
from utils.graphing_tools import graph_rate

group_column = "sample_group"

class RateCalculator():
    def __init__(self, model_path, out_path, graph_folder, settings_path, biomolecule_type):
        settings.load(settings_path)
        if model_path[-4:] == ".tsv":
            self.model = pd.read_csv(model_path, sep='\t')
        elif model_path[-4:] == ".csv":
            self.model = pd.read_csv(model_path)
        else: #$should never trigger unless we are fiddling with the gui
            raise ValueError("invalid file extension")
        self.out_path = out_path
        self.rate_model = None
        self.datapoint_model = None  # CQ
        self.graph_folder = graph_folder
        self.settings_path = settings_path
        #$get the number of cores we're using for multiprocessing
        if settings.recognize_available_cores is True:
            self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
        #$breaks windows/python interactions if too many cores are used. very niche application but still relevant
        if self._n_processors > 60:
            self.n_processors = 60
        self._mp_pool = mp.Pool(self._n_processors)

        self.biomolecule_type = biomolecule_type  # CQ

    def write(self):
        self.rate_model.to_csv(
            path_or_buf=self.out_path,
            index=False
        )
        self.datapoint_model.to_csv(  # CQ
            path_or_buf=self.out_path[:-4] + "_datapoints.tsv",
            sep='\t',
            index=False
        )
    
    #$this function governs which types of calculations to do and which
    #$ rate equation to use. collects results and puts them in self.rate_model
    def calculate(self):
        #$ catch the optimize warnings and so on
        w.filterwarnings("error")
        rate_results = []
        datapoint_results = []  # CQ
        max_time = max(self.model["time"])
        if self.biomolecule_type == "Peptide":
            id_col = settings.peptide_analyte_id_column
        elif self.biomolecule_type == "Lipid":
            id_col = settings.lipid_analyte_id_column
        else:
            print("Unknown Analyte Type")
            return
        #$could put manual bias correction here, would be slightly more 
        #$efficient, but make the logic less clear
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
            
        #$we're going to load up a function for mp.  we need to load up the 
        #$function. don't add fn_col, fn_std_dev, calc_type, and manual_bias,
        #$they will be added later
        rate_function = partial(RateCalculator._mp_function, 
                                settings_path = self.settings_path,
                                graph_folder = self.graph_folder,
                                rate_eq = rate_eq,
                                max_time = max_time,
                                p0 = p0)
            
        #$call the rate for each relevant measurement type
        #$add measurement specific arguments to partial function
        if settings.use_abundance != "No":
            temp_rate_function = partial(rate_function,
                                         calc_type = "Abundance",
                                         fn_col = "abund_fn",
                                         fn_std_dev = "frac_new_abunds_std_dev",
                                         manual_bias = settings.abundance_manual_bias,
                                         std_dev_filter = settings.abundance_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                results = list(
                    tqdm(
                        self._mp_pool.imap_unordered(temp_rate_function, groups),
                        total=len(groups), desc="Abundance Rate Calculation: "
                    )
                )
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Abundance Rate Calculation: "):
                    results.append(temp_rate_function(group))
            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ
            
        if settings.use_neutromer_spacing:
            temp_rate_function = partial(rate_function,
                                         calc_type = "Spacing",
                                         fn_col = "nsfn",
                                         fn_std_dev = "frac_new_mzs_std_dev",
                                         manual_bias = settings.spacing_manual_bias,
                                         std_dev_filter = settings.spacing_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                results = list(
                    tqdm(
                        self._mp_pool.imap_unordered(temp_rate_function, groups),
                        total=len(groups), desc="Spacing Rate Calculation: "
                    )
                )
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Spacing Rate Calculation: "):
                    results.append(temp_rate_function(group))
            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ

        if settings.use_abundance != "No" and settings.use_neutromer_spacing:
            temp_rate_function = partial(rate_function,
                                         calc_type = "Combined",
                                         fn_col = "cfn",
                                         fn_std_dev = "frac_new_combined_std_dev",
                                         manual_bias = settings.combined_manual_bias,
                                         std_dev_filter = settings.combined_agreement_filter,
                                         biomolecule_type=self.biomolecule_type)
            results = list()
            if settings.debug_level == 0:
                results = list(
                    tqdm(
                        self._mp_pool.imap_unordered(temp_rate_function, groups),
                        total=len(groups), desc="Combined Rate Calculation: "
                    )
                )
            elif settings.debug_level >= 1:
                for group in tqdm(groups, desc="Combined Rate Calculation: "):
                    results.append(temp_rate_function(group))
            rate_results.append(pd.DataFrame([i[0] for i in results]))  # CQ
            datapoint_results.append([i[1] for i in results])  # CQ
            
 
        #$ from https://stackoverflow.com/questions/44327999/python
        #$-pandas-merge-multiple-dataframes answer 1 accessed 9/16/2020
        #$merges any number of dfs into one df
        self.rate_model = reduce(lambda  left,right: pd.merge(left,right,
                                on=['analyte_id', 'analyte_name', "group_name"],
                                how='outer'), rate_results)
        self.datapoint_model = pd.concat(datapoint_results[0])
        #$swap back to normal warning behaviour
        w.filterwarnings("default")
        
        self._mp_pool.close()
        self._mp_pool.join()
            
        if not settings.verbose_rate:
            self.trim_verbose_data()
        
    def trim_verbose_data(self):
        needed_columns = ['analyte_id', 'analyte_name', "group_name"]
        potential_columns = ["{} rate", "{} 95pct_confidence", "{} half life"]
        if settings.asymptote != "Fixed":
            potential_columns.insert(1, "{} asymptote")
        if settings.use_abundance:
            for p in potential_columns:
                needed_columns.append(p.format("Abundance"))
        if settings.use_neutromer_spacing:
            for p in potential_columns:
                needed_columns.append(p.format("Spacing"))
        if settings.use_abundance and settings.use_neutromer_spacing:
            for p in potential_columns:
                needed_columns.append(p.format("Combined"))
        self.rate_model = self.rate_model[needed_columns]

    def _mp_function(data_tuple, settings_path, fn_col, 
                     fn_std_dev, calc_type, manual_bias, std_dev_filter, graph_folder,
                     rate_eq, max_time, p0, biomolecule_type):
        w.filterwarnings("error")
        settings.load(settings_path)
        pd.options.mode.chained_assignment = None
        
        id_values, group = data_tuple[0], data_tuple[1]
        id_name = id_values[0]
        sample_group_name =id_values[1]
        if biomolecule_type == "Peptide":
            common_name = group[settings.peptide_analyte_name_column].iloc[0]
        elif biomolecule_type == "Lipid":
            common_name = group[settings.lipid_analyte_name_column].iloc[0]
        else:
            raise Exception("Not a Valid Biomolecule Type")
        #$drop error string could do earlier for more speed, but this is
        #$clearer and allows errors that affect only one calculation type
        group = RateCalculator._error_trimmer(
            group, [fn_col, fn_std_dev])
        #$the copy is just to avoid a SettingWithCopy warning in a few
        #$operations.  if it causes problems remove and suppress warning
        group = group[group[fn_std_dev] < std_dev_filter].copy()
        
        if len(group) == 0:
            result = RateCalculator._make_error_message(
                "No Isotope Envelopes Agree","", id_name, common_name,
                sample_group_name, calc_type, 0, 0, 0, 0)
            return result, group
    
        #offset all values by a certain amount (instrument bias)
        if settings.bias_calculation == "calculated":
            bias = RateCalculator._calc_bias(group, fn_col)
            group[fn_col] = group[fn_col] - bias
        elif settings.bias_calculation == "manual": #$ user designated bias
            group[fn_col] = group[fn_col] - manual_bias
        
        xs = np.concatenate(([0], group['time'].to_numpy()))

        ys = np.concatenate(([settings.y_intercept_of_fit],
                                 group[fn_col].to_numpy()))
        
            
        if settings.roll_up_rate_calc:
                xs, ys, devs = RateCalculator._roll(xs, ys)
        else:
            devs = np.concatenate((
                [settings.error_of_zero],
                group[fn_std_dev].to_numpy()))
            
        # Get the number of unique time points, and continue if not enough
        num_unique_times = len(set(group['time']))

        if biomolecule_type == "Peptide":
            unique_length = len(set(group[settings.unique_sequence_column]))
        # $for lipids or metaboloites or similar, the unique length means
        # $nothing.  if needed can add the elif
        else:
            unique_length = ""
        num_measurements = len(group.index)
        
        #TODO$ this is not ideal but this is a good first attempt
        num_files = len(set(group["mzml_path"]))
        
        #TODO: Handle Num Bio Reps in Stuff
        
        # I think this fixes the issues with technical replicates.
        #$need to use astype(str) on all or if someone uses numbers for group names or replicate names issues result
        num_bio_reps = len(set(group["time"].astype(str) + group["sample_group"].astype(str) + group["bio_rep"].astype(str)))
        
        if num_unique_times < settings.minimum_nonzero_points:
            result = RateCalculator._make_error_message(
                "Insufficient Timepoints","", id_name, common_name, sample_group_name,
                calc_type, num_measurements, num_unique_times, unique_length,
                num_files)
            return result, group
        # perform fit
        try:
            #$DO NOT use std dev as the Sigma because it creates influential outliers
            #$don't use sigma unless we have a different
            popt, pcov = curve_fit(
                f=rate_eq, xdata=xs, ydata=ys,
                p0=p0)
            
             # pull results of fit into variables
            rate = popt[0]
            asymptote = \
                popt[1] if len(popt) > 1 else settings.fixed_asymptote_value
            #TODO$ ci uses degrees of freedom = n-k where n is the number of points and k is the number of parameters estimated
            #$including intercept in linear regression.  if asymptote is fixed k =1 otherwise k =2 (intercept is fit by equation, not data)
            #$not counting charge states and different peptides as unique measurements.
            #$despite the claim in the documentation, acorrding to statistics consultation and every site I checked, the np.sqrt(np.diag(pcov))[0]
            #$is standard error, not std dev, so don't divide by sqrt of n
            
            confint = \
                t.ppf(.975, num_files - len(popt)) * \
                np.sqrt(np.diag(pcov))[0]
                
            y_predicted = dur.simple(xs, rate, asymptote, settings.proliferation_adjustment)
            r_2 = dur.calculate_r2(ys, y_predicted)
            
            result = {
                'analyte_id': id_name,
                'analyte_name': common_name,
                'group_name': sample_group_name,
                '{} rate'.format(calc_type) : rate,
                '{} asymptote'.format(calc_type) : asymptote,
                '{} std_error'.format(calc_type): np.sqrt(np.diag(pcov))[0],
                '{} 95pct_confidence'.format(calc_type): confint,
                '{} half life'.format(calc_type): RateCalculator._halflife(rate),
                '{} R2'.format(calc_type): r_2,
                "{} files observed in".format(calc_type): num_files,
                '{} num_measurements'.format(calc_type):
                    num_measurements,
                '{} num_time_points'.format(calc_type):
                    num_unique_times,
                '{} uniques'.format(calc_type): unique_length,
                '{} exceptions'.format(calc_type): "",
                #$'calculation_type': calc_type
            }
            #$ if there is an asymptote need to provide it
            if biomolecule_type == "Peptide":
                graph_name = "{}_{}_{}".format(id_name, sample_group_name,
                                               fn_col)
                graph_title = "{}_{}_{}\nk={}, a={}".format(common_name,
                                                            sample_group_name, fn_col, result[f'{calc_type} rate'],
                                                            1.0)
               if settings.roll_up_rate_calc:
                  graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                           rate_eq, graph_folder, max_time,
                           settings.asymptote, devs, biomolecule_type,title=graph_title)
               else:
                  graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                           rate_eq, graph_folder, max_time,
                           settings.asymptote, biomolecule_type, title=graph_title)
            elif biomolecule_type == "Lipid":
                graph_name = "{}_{}_{}".format(common_name,
                                               sample_group_name, fn_col)
                graph_title = "{}_{}_{}\nk={}, a={}".format(common_name,
                                                            sample_group_name, fn_col, result['Abundance rate'],
                                                            result['Abundance asymptote'])
               if settings.roll_up_rate_calc:
                  graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                           rate_eq, graph_folder, max_time,
                           settings.asymptote, devs, biomolecule_type, full_data=group, title=graph_title)
               else:
                  graph_rate(graph_name, xs, ys, rate, asymptote, confint,
                           rate_eq, graph_folder, max_time,
                           settings.asymptote, biomolecule_type, full_data=group, title=graph_title)
            else:
                graph_name = None
                graph_title = None
            
            
        except Exception as c:
            #$"we have a guess but are unsure" warning
            if type(c).__name__ == "OptimizeWarning":
                current_exception = \
                    'OptimizeWarning: optimal fit could not be found'
            #$couldn't find the minimum
            elif type(c).__name__ == "RuntimeError":
                current_exception = \
                    'fit could not be found'
            else:
                raise c #$will stop here so don't need to consider further
            result = RateCalculator._make_error_message(
                "value could not be determined", current_exception, id_name,
                common_name, sample_group_name, calc_type,num_measurements,
                num_unique_times, unique_length, num_files)
        
        return result, group
            
    #$funciton that does a mad oulier check on a numpy array
    @staticmethod
    def _mask_outliers(y_values):
        if len(y_values == 1):return y_values
        med = np.median(y_values)
        differences = abs(y_values - med)
        med_abs_dev = np.median(differences)
        z_values = (differences * .6745) / med_abs_dev
        good_indicies = np.where(z_values < settings.zscore_cutoff)[0]
        return y_values[good_indicies]

    #$roll up time points 
    @staticmethod
    def _roll(old_x, old_y):
        new_x = np.unique(old_x)
        new_y, errors = [], []
        for x in new_x:
            #$ grab indicies for the time
            x_time_indicies = np.where(old_x == x)[0]
            #$grab only relevant indicies of y and oulier check
            temp_y_values = RateCalculator._mask_outliers(
                old_y[x_time_indicies])
            #$second length check (first is in _mask_outliers) to avoid
            #$medians and std devs on statistically invalid data
            if len(temp_y_values) == 1:
                new_y.extend(temp_y_values)
                errors.append(settings.error_of_non_replicated_point)
            else:
                new_y.append(np.median(temp_y_values))
                standard_dev= np.std(temp_y_values, ddof =1)
                #$I have seen std dev be 0 which causes problems for the fit, 
                #$so we'll just force no 0s
                if standard_dev == 0.0: errors.append(settings.error_of_zero)
                else: errors.append(standard_dev)
        return new_x, np.array(new_y), np.array(errors)
            
    
    #$need to calculate median value at time zero in the target column
    #$already checked for std dev error at this point
    @staticmethod
    def _calc_bias(df, target_column):
        test_df = df[df["time"] == 0]
        if len(test_df) > 0:
            return np.median(test_df[target_column])
        else:
            return 0
    
    #$implemented and working
    @staticmethod
    def _halflife(x):
        return np.log(2)/x
    
    #$just to shorten the code for the error messages
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
                  '{} uniques'.format(calc_type):unique_length,
                  '{} exceptions'.format(calc_type): exception
                  }
        return result
    
    #$removes errors from a previous analysis step, which are strings
    #$just cut out those rows
    #$ df is the dataframe to trim, error columns is a list of columns where
    #$errors are located. returns the trimmed dataframe
    @staticmethod
    def _error_trimmer(df, error_columns):
        df[error_columns] = df[error_columns].apply(
                pd.to_numeric, errors = 'coerce')
        return df.dropna(axis =0, subset = error_columns)
