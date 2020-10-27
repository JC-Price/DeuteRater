'''Calculation of turnover rates

'''

from functools import partial, reduce
import pandas as pd
import numpy as np
from scipy.stats import t  # pearsonr?
from scipy.optimize import curve_fit  # least_sq?
import warnings as w

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.rate_equations as dur
from utils.graphing_tools import graph_rate

group_column = "sample_group"

class RateCalculator():
    def __init__(self, model_path, out_path, graph_folder, settings_path,
                 biomolecule_type):
        settings.load(settings_path)
        self.model = pd.read_csv(model_path, sep ="\t")
        self.out_path = out_path
        self.rate_model = None
        self.graph_folder = graph_folder
        self.biomolecule_type = biomolecule_type

    def write(self):
        self.rate_model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )
    
    #$this function governs which types of calculations to do and which
    #$ rate equation to use. collects results and puts them in self.rate_model
    def calculate(self):
        #$ catch the optimize warnings and so on
        w.filterwarnings("error")
        rate_results = []
        max_time = max(self.model["time"])
        if self.biomolecule_type == "Peptide":
            id_col = settings.peptide_analyte_id_column
        elif self.biomolecule_type == "Lipid":
            id_col = settings.lipid_analyte_id_column
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
        #$call the rate for each relevant measurement type
        if settings.use_abundance:
            rate_results.append(self.rate_for_calc_type(
                groups, p0, "Abundance", "afn", "frac_new_abunds_std_dev", 
                rate_eq, settings.abundance_agreement_filter,
                settings.abundance_manual_bias, max_time))
        if settings.use_neutromer_spacing:
            rate_results.append(self.rate_for_calc_type(
                groups, p0, "Spacing", "nsfn", "frac_new_mzs_std_dev", 
                rate_eq, settings.spacing_agreement_filter,
                settings.spacing_manual_bias, max_time))
        if settings.use_abundance and settings.use_neutromer_spacing:
            rate_results.append(self.rate_for_calc_type(
                groups, p0, "Combined", "cfn", "frac_new_combined_std_dev", 
                rate_eq, settings.combined_agreement_filter,
                settings.combined_manual_bias, max_time))
        #$ from https://stackoverflow.com/questions/44327999/python
        #$-pandas-merge-multiple-dataframes answer 1 accessed 9/16/2020
        #$merges any number of dfs into one df
        self.rate_model = reduce(lambda  left,right: pd.merge(left,right,
                                on=['analyte_id', 'analyte_name'], 
                                how='outer'), rate_results)
        #$swap back to normal warning behaviour
        w.filterwarnings("default")
        
    #$actually does the rate calculation and returns results or an error
    def rate_for_calc_type(self, groups,p0, calc_type, fn_col, fn_std_dev, 
                           rate_eq, std_dev_filter, manual_bias, max_time):
        calc_results = []
        # fit each group
        for id, group in tqdm(groups):
            id_name = id[0]
            sample_group_name =id[1]
            if self.biomolecule_type == "Peptide":
                common_name = group[settings.peptide_analyte_name_column].iloc[0]
            if self.biomolecule_type == "Lipid":
                common_name = group[settings.lipid_analyte_name_column].iloc[0]
            #$drop error string could do earlier for more speed, but this is 
            #$clearer and allows errors that affect only one calculation type
            group = RateCalculator._error_trimmer(
                group, [fn_col, fn_std_dev])
            #$the copy is just to avoid a SettingWithCopy warning in a few
            #$operations.  if it causes problems remove and suppress warning
            group = group[group[fn_std_dev] < std_dev_filter].copy()
            if len(group) == 0:
                result = self._make_error_message(
                    "No Isotope Envelopes Agree","", id_name, common_name,  
                    sample_group_name, calc_type, 0, 0, 0, 0)
                calc_results.append(result)
                continue
            
            #offset all values by a certain amount (instrument bias)
            if settings.bias_calculation == "calculated":
                bias = self._calc_bias(group, fn_col)
                group[fn_col] = group[fn_col] - bias
            elif settings.bias_calculation == "manual": #$ user designated bias  
                group[fn_col] = group[fn_col] - manual_bias
                
            xs = np.concatenate(([0], group['time'].to_numpy()))

            ys = np.concatenate(([settings.y_intercept_of_fit], 
                                 group[fn_col].to_numpy()))
            
            if settings.roll_up_rate_calc:
                xs, ys, devs = self._roll(xs, ys)
            else:
                devs = np.concatenate((
                    [settings.error_of_zero],
                    group[fn_std_dev].to_numpy()))
            
            # Get the number of unique time points, and continue if not enough
            num_unique_times = len(set(group['time']))
            if self.biomolecule_type == "Peptide":
                unique_length = len(set(group[settings.unique_sequence_column]))
            #$for lipids or metaboloites or similar, the unique length means
            #$nothing.  if needed can add the elif
            else:
                unique_length = ""
            num_measurements = len(group.index)
            #TODO$ this is not ideal but this is a good first attempt
            num_files = len(set(group["mzml_path"]))
            
            if num_unique_times < settings.minimum_nonzero_points:
                result = self._make_error_message(
                    "Insufficient Timepoints","", id_name, common_name, sample_group_name,
                    calc_type, num_measurements, num_unique_times, unique_length,
                    num_files)
                
                calc_results.append(result)
                continue

            # perform fit
            try: 
                popt, pcov = curve_fit(
                    f=rate_eq, xdata=xs, ydata=ys,
                    p0=p0, sigma=devs
                )
    
                # pull results of fit into variables
                rate = popt[0]
                asymptote = \
                    popt[1] if len(popt) > 1 else settings.fixed_asymptote_value
                #TODO$ may need to adjust the ci value to calculate, but for 
                #$now num_files works
                confint = \
                    t.ppf(.975, num_files - 1) * \
                    np.sqrt(np.diag(pcov))[0] / \
                    np.sqrt(num_files)
                result = {
                    'analyte_id': id_name,
                    'analyte_name': common_name,
                    'group_name': sample_group_name,
                    '{} rate'.format(calc_type) : rate,
                    '{} asymptote'.format(calc_type) : asymptote,
                    '{} std_dev'.format(calc_type): np.sqrt(np.diag(pcov))[0],
                    '{} 95pct_confidence'.format(calc_type): confint,
                    '{} half life'.format(calc_type): self._halflife(rate),
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
                if self.biomolecule_type == "Peptide":
                    graph_name = "{}_{}_{}".format(id_name, sample_group_name, 
                                                   fn_col)
                elif self.biomolecule_type == "Lipid":
                    graph_name = "{}_{}_{}".format(common_name, 
                                                   sample_group_name, fn_col)
                if settings.roll_up_rate_calc:
                    graph_rate(graph_name, xs, ys, rate, asymptote, confint, 
                               rate_eq, self.graph_folder, max_time, 
                               settings.asymptote, devs)
                else:
                    graph_rate(graph_name,  xs, ys, rate, asymptote, confint, 
                               rate_eq, self.graph_folder, max_time, 
                               settings.asymptote)
                               
                               
            
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
                result = self._make_error_message(
                    "value could not be determined", current_exception, id_name, 
                    common_name, sample_group_name, calc_type,num_measurements,  
                    num_unique_times, unique_length, num_files)
            
            calc_results.append(result)
           
        return pd.DataFrame(calc_results)

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
                  '{} std_dev'.format(calc_type): main_error,
                  '{} 95pct_confidence'.format(calc_type): main_error,
                  '{} half life'.format(calc_type): main_error,
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
    """
    #$function has been moved to the fraction new creator
    #$if need be we can turn that off and use this instead
    #$ this function has not been tested so would need some work to 
    def _summarize_acfn(self, calculation_type):
        # TODO: this needs to add columns, not just be a simple eq
        fraction_news = []
        # TODO: make sure to deal with outliers
        if (calculation_type == 'abundance'):
            fraction_news.extend(self.model['frac_new_abunds'].tolist())
        elif calculation_type == 'spacing':
            fraction_news.extend(self.model['frac_new_mzs'].tolist())
        elif calculation_type == 'combined':
            fraction_news.extend(
                self.model['frac_new_abunds'].tolist() +
                self.model['frac_new_mzs'].tolist()
            )
        summaries = [mean, stdev, len]
        return [fn(fraction_news) for fn in summaries]
    """