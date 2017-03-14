"""
Copyright (c) 2016 Bradley Naylor, Michael Porter, J.C. Price, and Brigham Young University

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
import numpy as np
import pandas as pd
from scipy.stats import t
from scipy.optimize import curve_fit
import os

def remove_string(df, col, string, all_cols, numeric = True):
    if df[col].dtype != np.float64:
        new_df = df[df[col] != string].copy()
        for column in all_cols:
            if numeric:
                new_df.loc[:,column] = pd.to_numeric(new_df[column])
        return new_df
    else:
        return df


def remove_outlier(turnovers, filters):
    med = np.median(turnovers)
    diff = []
    for point in  turnovers:
        diff.append(abs(point-med))
    mad = np.median(diff)
    zScore = []
    for item in diff:
        z = .6745 *item/mad
        zScore.append(z)
    true_turnovers =[]
    positions = []
    for pos in range(len(zScore)):
        if zScore[pos] < filters["ZCUTOFF"]:
            true_turnovers.append(turnovers[pos])
            positions.append(pos)
    return true_turnovers, positions

#compresses data and filters on stdev and rt_diff.  if filtering needs to be changed, do it here
def mida_compressor(data, relevant_change, measured, filters):
    pre_compress = data.groupby(['time', 'Protein ID'])
    compressed_data= []
    for name, group in pre_compress:
        protein_name = list(group["Protein name"])[0]
        if len(group[relevant_change]) >1:
            before_outlier = len(group[relevant_change])
            new_group, positions = remove_outlier(list(group[relevant_change]), filters)
            after_outlier = len(new_group)
            sequences = []
            for p in positions:
                sequences.append(list(group['sequence'])[p])
            if len(new_group) ==1: x = 0 #std. dev occasionally gives values for things of length 1 (rounding error) This forces that not to happen
            else: x = np.std(new_group, ddof =1)
            compressed_data.append([name[0], name[1], protein_name, np.median(new_group),x , len(new_group), len(set(sequences)),  before_outlier, after_outlier, list(set(sequences))])
        else:
            if len(group) ==1: x = 0 #std. dev occasionally gives values for things of length 1 (rounding error) This forces that not to happen
            else: x = np.std(group[relevant_change], ddof =1)
            unique_peptides = group['sequence'].unique()
            compressed_data.append([name[0], name[1], protein_name, np.median(group[relevant_change]),x , len(group), len(unique_peptides),"only one measurement. no outlier check","only one measurement. no outlier check",  unique_peptides])
    return pd.DataFrame(compressed_data, columns = ['time', 'Protein ID', "Protein name", '{0} fraction new'.format(measured), '{0} std_dev'.format(measured), 'number of measurements', 'number of unique peptides', "# of peptides before outlier removal", "# of peptides after outlier removal", 'unique peptides' ])

#does the actual rate calculation
def mida_rates(data,water, measurement, filters, settings):
    import warnings as w
    w.filterwarnings("error")
    #the water dataframe has the time as well as the water so this adds the time (in the other function we only used water percentage
    if not "time" in data.columns:
        def add_time(row):
            row['time'] = water['time'][row['file']]
            return row
        data = data.apply(add_time, axis =1)
    data = data[data['time'] != 0] #filter step to remove the zero hour timepoint.  since there is no labeling at time zero(and the fit equation will not allow it) this would only add noise
    pre_compress = data.groupby('Protein ID')
    rate_data = []
    if settings["asymptote"] == "Variable":
        std_dev_loc = 1
        def rate_equation(t,a, k):
            return a - a*(np.exp(-(k+settings["Proliferation Adjustment"])*t))
        p0 = [1,.1]
    else:
        std_dev_loc = 0
        def rate_equation(t,k):
            return settings["fixed"] - settings["fixed"]*(np.exp(-(k+settings["Proliferation Adjustment"])*t))
        p0 = .1
    #these choose the peak to use.  alter for the appropriate header if change is needed
    if settings["Roll Up"]:
        y_loc = '{0} fraction new'.format(measurement)
    elif measurement == "abundance":
        y_loc = 'M0 fraction new for calculation'
    else:
        y_loc ='median {0} fraction new'.format(measurement)
    for name, group in pre_compress:
        if len(set(group['time'])) >= settings["Minimum Non-Zero Points"]:
            protein_name = list(group["Protein name"])[0]
            #set zeroes
            x_values = [0]
            if settings["Roll Up"]:
                error_values = [filters["Error of Zero"]]
            y_values = [0]
            #add actual values to the x, y and errors(if the roll_up is present)
            x_values.extend(group['time'])
            y_values.extend(group[y_loc])
            if settings["Roll Up"]:
                error_values.extend(group['{0} std_dev'.format(measurement)])
                for e in range(len(error_values)):
                    #error for only one point should be assumed to be the same as for zero (instrument background)
                    if error_values[e] == 0:
                        error_values[e] = filters["Error of non-replicated point"]
                unique_peptides = []
                for data in group['unique peptides']:
                    unique_peptides.extend(data)
                try:
                    k, pcov = curve_fit(rate_equation, np.asarray(x_values), np.asarray(y_values), p0, np.asarray(error_values))
                    if settings["asymptote"] == "Variable": 
                        asymptote = k[0]
                        rate = k[1]
                    else: 
                        asymptote = settings["fixed"]
                        rate = k[0]
                    final_row = [name, protein_name, rate, asymptote, np.sqrt(np.diag(pcov))[std_dev_loc], t.ppf(.975, sum(group['number of measurements'])-1)* np.sqrt(np.diag(pcov))[std_dev_loc]/np.sqrt(sum(group['number of measurements'])), sum(group['number of measurements']), len(set(unique_peptides)), '']
                    if settings["asymptote"] == "Variable": final_row.append(np.sqrt(np.diag(pcov))[0]) 
                    rate_data.append(final_row)
                #way of unpacking excpetion name from http://stackoverflow.com/questions/9823936/python-how-do-i-know-what-type-of-exception-occured answer 1 accessed 2/11/2016
                except Exception as c:
                    if type(c).__name__ == "OptimizeWarning":
                        final_row = [name, protein_name, "k could not be determined", "a could not be determined", "Std. Dev. of k could not be determined", "95% C.I. could not be determined",sum(group['number of measurements']), len(set(unique_peptides)), 'OptimizeWarning: optimal fit could not be found']
                        if settings["asymptote"] == "Variable": final_row.append("std. dev. of a could not be determined") 
                        rate_data.append(final_row)
                    elif type(c).__name__ == "RuntimeError":
                        final_row = [name, protein_name, "k could not be determined", "a could not be determined", "Std. Dev. of k could not be determined", "95% C.I. could not be determined",sum(group['number of measurements']), len(set(unique_peptides)), 'fit could not be found']
                        if settings["asymptote"] == "Variable": final_row.append("std. dev. of a could not be determined") 
                        rate_data.append(final_row)
                    else: 
                        raise c
            else:
                measurements = len(group['sequence'])
                unique_peptides = len(set(group['sequence']))
                try:
                    k, pcov = curve_fit(rate_equation, np.asarray(x_values), np.asarray(y_values), p0)
                    if settings["asymptote"] == "Variable": 
                        asymptote = k[0]
                        rate = k[1]
                    else: 
                        asymptote = settings["fixed"]
                        rate = k[0]
                    #adjust for t with the std error for 95% C.I.
                    final_row = [name, protein_name, rate, asymptote, np.sqrt(np.diag(pcov))[std_dev_loc], t.ppf(.975, measurements-1)* np.sqrt(np.diag(pcov))[std_dev_loc]/np.sqrt(measurements), measurements, unique_peptides, '']
                    if settings["asymptote"] == "Variable": final_row.append(np.sqrt(np.diag(pcov))[0])
                    rate_data.append(final_row)
                except Exception as c:
                    if type(c).__name__ == "OptimizeWarning":
                        final_row = [name, protein_name, "k could not be determined", "a could not be determined", "Std. Dev. of k could not be determined", "95% C.I. could not be determined", measurements, unique_peptides, 'OptimizeWarning: optimal fit could not be found']
                        if settings["asymptote"] == "Variable": final_row.append("Std. Dev. of a could not be determined")
                        rate_data.append(final_row)
                    elif type(c).__name__ == "RuntimeError":
                        final_row = [name, protein_name, "k could not be determined", "a could not be determined", "Std. Dev. could not be determined", "95% C.I. could not be determined", measurements, unique_peptides, 'a fit could not be found'] 
                        if settings["asymptote"] == "Variable": final_row.append("Std. Dev. of a could not be determined")
                        rate_data.append(final_row)
                    else: 
                        raise c
    #return all the data if successful, otherwise returns a string to tell the main function there is an error (which subsequently tells the user)
    w.filterwarnings("default")
    if len(rate_data) > 0:
        header = ['Protein ID', "Protein name", '{0} rate'.format(measurement), "Asymptote",'{0} std. dev'.format(measurement), '{0} 95% confidence'.format(measurement), 'number of measurements', 'unique_peptides', "Problems"]
        if settings["asymptote"] == "Variable":
            header.append("a std. dev")
        return pd.DataFrame(rate_data, columns = header)
    else: return 'error'
def halflife(x):
    if type(x) == str:
        return "Half-life could not be determined"
    else: return np.log(2)/x

def rate_prep(data, analysis, settings, filters, water, Messages, Graph):
    lower_name = analysis.lower()
    if analysis == "Abundance": use_col = "M0 fraction new for calculation"
    else: use_col = "median {0} fraction new".format(lower_name)
    #remove strings and np.nan's and then filter (filter does not work with strings)
    data = data[data['{0} std_dev'.format(lower_name)].apply(np.isreal)]
    data_filtered = data[data['{0} std_dev'.format(lower_name)] <= filters["{0} Agreement Filter".format(analysis)]]
    if len(data_filtered > 0):
        #deal with roll_up or not for rate calculatiosn
        if settings["Roll Up"]:
            if "time" not in data_filtered.columns:
                def add_time(row):
                    row['time'] = water['time'][row['file']]
                    return row
                data_filtered = data_filtered.apply(add_time, axis =1)
            data_compressed = mida_compressor(data_filtered, use_col, lower_name, filters)
            data_compressed.loc[:,:"# of peptides after outlier removal"].to_csv(os.path.join(settings["Output Folder"], '{}_Protein_Roll_Up.csv'.format(analysis)), index = False)
            data_rates = mida_rates(data_compressed, water, lower_name, filters, settings)
        else:
            data_rates = mida_rates(data_filtered, water, lower_name,filters, settings)
        #make sure data is there, calculate half-life, and then print it out
        if type(data_rates) != str:
            data_rates['{0} half-life'.format(lower_name)] = data_rates['{0} rate'.format(lower_name)].apply(halflife)
            if settings["Roll Up"]:
                data_rates.to_csv(os.path.join(settings["Output Folder"], 'Final_{}_Rates_Roll_Up.csv'.format(analysis)), index = False)
            else:
                data_rates.to_csv(os.path.join(settings["Output Folder"], 'Final_{}_Rates.csv'.format(analysis)), index = False)
            Messages.append("{0} Analysis Completed Successfully".format(analysis))
            Graph[analysis] = True
        else:
            Messages.append("{0} data was not present at enough times to fit a rate for any protein. {0} Calculations will be ignored".format(analysis))
            Graph[analysis] = False
    else: 
        Messages.append("{0} data was not present at enough times to fit a rate for any protein. {0} Calculations will be ignored".format(analysis))
        Graph[analysis]= False
    return Messages, Graph

def calculate_bias(complete, calc_type, filters):
    if calc_type == "abundance":
        key_column = "M0 fraction new for calculation"
        cutoff = filters["Abundance Agreement Filter"]
    else:
        key_column = "median {} fraction new".format(calc_type)
        if calc_type == "combined": cutoff = filters["Combined Agreement Filter"]
        else: cutoff = filters["neutromer spacing Agreement Filter"]
    determination = complete[(complete["{} std_dev".format(calc_type)]  < cutoff) & (complete['Experimental_Amount_of_Label (eaol) (0 = {0})'.format(filters["Zero Labeling Check"])]==0)]
    if len(determination) > 0:
        return np.median(determination[key_column])
    else:
        return "no zero labeling for {} bias calculation. add data with zero labeling, relax std dev filters or change bias settings".format(calc_type)
def do_rates(complete_with_zero, complete_no_zero, calc_type, settings, filters, Messages, Graph, water):
    pd.options.mode.chained_assignment = None
    if calc_type == "abundance": 
        key_column =  "M0 fraction new for calculation"
        complete_no_zero = remove_string(complete_no_zero, "{} std_dev".format(calc_type), "M0 could not be used", ["{} std_dev".format(calc_type), key_column], False)
    else: key_column = "median {} fraction new".format(calc_type)
    if calc_type == "neutromer spacing": filter_name = "neutromer spacing"
    else: filter_name = calc_type.capitalize()
    continue_calc = True
    improved_data = remove_string(complete_no_zero, "{} std_dev".format(calc_type), "insufficient labeling", ["{} std_dev".format(calc_type), key_column])
    if settings["bias_correction"]=="DeuteRater Calculated":
        bias = calculate_bias(remove_string(complete_with_zero,"{} std_dev".format(calc_type), "insufficient labeling", ["{} std_dev".format(calc_type), key_column]) , calc_type, filters)
        if type(bias) ==str:
            continue_calc = False
            Messages.append(bias)
        else:
            improved_data[key_column] = improved_data[key_column]-bias
            improved_data.to_csv(os.path.join(settings["Output Folder"], '{}_bias_corrected.csv'.format(calc_type)), index = False, columns = ["file", "Protein ID", "Protein name", "sequence", key_column,"{} std_dev".format(calc_type)])
    if settings["bias_correction"] == "User Defined":
        improved_data[key_column] = improved_data[key_column]-settings["{}_bias".format(calc_type)]
        improved_data.to_csv(os.path.join(settings["Output Folder"], '{}_bias_corrected.csv'.format(calc_type)), index = False, columns = ["file", "Protein ID", "Protein name", "sequence", key_column,"{} std_dev".format(calc_type)])
    if continue_calc:
        Messages, Graph = rate_prep(improved_data, filter_name, settings, filters, water, Messages, Graph)
    return Messages, Graph

def rate_calc(Complete, water, settings, filters):
    Messages = []
    Graph = {"Abundance":False, "neutromer spacing": False, "Combined": False}
    Complete_no_zero = Complete[Complete['Experimental_Amount_of_Label (eaol) (0 = {0})'.format(filters["Zero Labeling Check"])]>0] # no point calculating with no label as it is a special case for calculation, and negative values make no sense
    if len(Complete_no_zero) > 0:
        #these values may need to change.  we don't want to change the main settings in case of further analysis
        if settings["Use Abundance"]:
            do_rates(Complete, Complete_no_zero, "abundance", settings, filters, Messages, Graph, water)
        if settings["Use neutromer spacing"]:
            do_rates(Complete, Complete_no_zero, "neutromer spacing", settings, filters, Messages, Graph, water)
        #$ last function aside from mida_rates
        if Graph["Abundance"] and Graph["neutromer spacing"]:
            do_rates(Complete, Complete_no_zero, "combined", settings, filters, Messages, Graph, water)
        elif Graph["Abundance"] != settings["Use Abundance"] or Graph["neutromer spacing"] != settings["Use neutromer spacing"]:
            Messages.append("Loss of Abundance or neutromer spacing due to insufficient data prevents combined analysis")
        if not (Graph["Abundance"] or Graph["neutromer spacing"]):
            Messages.append("No peaks agreed well enough for rate analysis.  Get better agreement or relax the filters")
    return Messages, Graph
