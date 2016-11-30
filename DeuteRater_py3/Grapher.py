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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import Rate_Calculator_vrs6 as rc


def bare_bones_graph(rate, asymptote,error, maximum, settings):
    plt.xlabel('Time')
    plt.ylabel('Fraction of Total Protein that is new')
    x = np.arange(0,maximum + maximum/10, .1)
    y = []
    y_plus_error = []
    y_minus_error = []
    def rate_equation(k,a, t):
        return a- a*np.exp(-(k+settings["Proliferation Adjustment"])*t)
    for i in x:
        y.append(rate_equation(rate, asymptote, i))
        plus = rate_equation(rate+error,asymptote, i)
        minus = rate_equation(rate-error,asymptote, i)
        if plus == float("inf") or plus > 10: plus = 10
        if minus == float("inf") or minus < 0: minus = 0
        y_plus_error.append(plus)
        y_minus_error.append(minus)
    plt.plot(x, y, 'k-')
    plt.plot(x, y_plus_error, 'k--')
    plt.plot(x, y_minus_error, 'k--')

def Make_Graph(analysis_type, input_folder, settings, filters):
    #$this is necessary due to inconsistent capitalization.  revise code if time permits
    lower_analysis = analysis_type.lower()
    if settings["Roll Up"]:
        data_file = os.path.join(input_folder,"{}_Protein_Roll_Up.csv".format(analysis_type))
        data_column = "{0} fraction new".format(lower_analysis)
        rate_file = 'Final_{}_Rates_Roll_Up.csv'.format(analysis_type)
    else:
        data_file = os.path.join(input_folder,"Calculation_of_Fraction_New_Protein.csv")
        rate_file = 'Final_{}_Rates.csv'.format(analysis_type)
        if analysis_type == "Abundance":
            data_column = "M0 fraction new for calculation"
        else:
            data_column = 'median {0} fraction new'.format(lower_analysis)
    output_folder = os.path.join(input_folder,"{0}_Graphs".format(analysis_type))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    """else: clear directory?"""
    time = pd.read_csv(os.path.join(input_folder,'Time_and_Labeling_Data.csv'))
    time = time.set_index(time["file"])
    max_time = time["time"].max()
    rates = pd.read_csv(os.path.join(input_folder,rate_file),usecols = ['Protein ID', '{0} rate'.format(lower_analysis), "Asymptote", '{0} 95% confidence'.format(lower_analysis)], low_memory = False)
    #http://stackoverflow.com/questions/22697773/how-to-check-the-dtype-of-a-column-in-python-pandas first answer (accessed 1/25/16)
    if rates["Asymptote"].dtype != np.float64:
        rates = rates[rates["Asymptote"] != "a could not be determined"]
        rates['{0} rate'.format(lower_analysis)] = rates['{0} rate'.format(lower_analysis)].astype(float)
        rates["Asymptote"] = rates["Asymptote"].astype(float)
        rates['{0} 95% confidence'.format(lower_analysis)] = rates['{0} 95% confidence'.format(lower_analysis)].astype(float)
    rates = rates.set_index(rates['Protein ID'])
    if settings["Roll Up"]:
        data = pd.read_csv(data_file, usecols = ['time', 'Protein ID', '{0} fraction new'.format(lower_analysis), '{0} std_dev'.format(lower_analysis)], low_memory = False)
    else:
        try:
            data = pd.read_csv(data_file, usecols = ['time', 'Protein ID', data_column, '{0} std_dev'.format(lower_analysis)], low_memory = False)
        except ValueError:
            data = pd.read_csv(data_file, usecols = ['file', 'Protein ID', data_column, '{0} std_dev'.format(lower_analysis)], low_memory = False)
            #these filters were applied previous to any roll up, so not needed above.  Are still needed here
            def add_time(row):
                row['time'] = time['time'][row['file']]
                return row
            data = data.apply(add_time, axis =1)
        if lower_analysis == "abundance": 
            filter_it_out = "insufficient labeling"
            # special error only in abund
            data = rc.remove_string(data, '{0} std_dev'.format(lower_analysis), "M0 could not be used", ['{0} std_dev'.format(lower_analysis), data_column], False)
        else: filter_it_out = "insufficient points"
        data = rc.remove_string(data, '{0} std_dev'.format(lower_analysis), filter_it_out, ['{0} std_dev'.format(lower_analysis), data_column])
        data = data[data['{0} std_dev'.format(lower_analysis)].apply(np.isreal)]
        data = data[data['{0} std_dev'.format(lower_analysis)] <= filters["{0} Agreement Filter".format(analysis_type)]]
    #$change to apply? if so index doesn't need to happen
    for i in rates["Protein ID"].tolist():
        bare_bones_graph(rates['{0} rate'.format(lower_analysis)][i], rates["Asymptote"][i], rates['{0} 95% confidence'.format(lower_analysis)][i], max_time, settings)
        points = data[data["Protein ID"] == i]#$ consider filtering out zeroes
        #rolled up proteins have error bars which are relevant.  without the roll up there are no error bars but far more points
        if settings["Roll Up"]:
            points = points.replace({'{0} std_dev'.format(lower_analysis): {0 :filters["Error of non-replicated point"]}})
            points.apply(lambda x: plt.errorbar(x["time"], x['{0} fraction new'.format(lower_analysis)], yerr = x['{0} std_dev'.format(lower_analysis)], elinewidth = 1, mfc = 'red', marker = 'o', linewidth = 0), axis =1)
        else:
            points.apply(lambda x: plt.plot(x["time"], x[data_column], 'ro'), axis = 1)
        plt.title(i)
        plt.savefig(os.path.join(output_folder,"{}.pdf".format(i)))
        plt.clf()