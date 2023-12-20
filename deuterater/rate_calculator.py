# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Christian Andersen, Chad Quilling, J.C. Price, and Brigham Young University
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

'''
This is the actual calculation of pepetide turnover rates and so is the largest
and most complex class. It is also the most expensive in terms of computer resources
and calculation time.  Be cautious running on a large dataset without multiprocessing
it will take a long time.

'''

import os
from functools import partial
import pandas as pd
import numpy as np
import scipy.integrate
from scipy.stats import binom
from scipy.optimize import brent
from scipy.integrate.odepack import ODEintWarning
from scipy.optimize import curve_fit
import scipy.interpolate as si
from collections import defaultdict



from pathlib import Path

import time
import warnings as w
import multiprocessing as mp
from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
from utils.graphing_tools import graph_rate_results, graph_optimization_of_error
from utils.emass import normalize

#$the spline fit used for the deuterium enrichment was previously used as a lambda with a fucntion call
#$but the pickler for multiprocessing will not take a lambda
#$this function does the spline calculations in a pickling safe manner.
class pickle_safe_spline(object):
    def __init__(self, time_values, deuterium_values):
        self.make_the_spline(time_values, deuterium_values)
    def make_the_spline(self, t, d):
        std_d = np.std(d)
        w_y = [3/std_d for i in range(len(t))]
        self.spl = si.splrep(t, d, w=w_y)
    #$adjusted to ensure it can't give negative values which can happen really early on
    #$ can't absolute value 
    def final_spline(self, x):
        return_value = si.splev(x, self.spl)
        if type(x) != np.ndarray:
            return max(return_value, np.asarray(0))
        else:
            return_value[return_value<0] =0
            return return_value

#$used for an error estimation graph
def parabola(x, a, b, c):
    return a*x**2 + b*x + c

group_column1 = "sample_id"
group_column2 = "Sequence"
time_column = "time"

#$need to deal with warnings.  in particular fitting warnings (ODEintWarning) should not be trusted
#$turn them into an error so we can catch.  will allow UserWarnings since it usually triggers on "dopri5: larger nsteps needed"
#$which we do not currently seek to catch. 
w.simplefilter("error", ODEintWarning) 
w.simplefilter("ignore", UserWarning) 


#$as with all the calculation steps this is a class for consistent calls in the main
class RateCalculator():
    def __init__(self, model_path, out_path, graph_folder_isotopes, graph_folder_optimization, settings_path):
       
        settings.load(settings_path)
        model_path = Path(model_path)
        #$allow .tsv and .csv
        if model_path.suffix == '.tsv':
            self.model = pd.read_csv(
                filepath_or_buffer=str(model_path),
                sep='\t'
                )
        elif model_path.suffix == '.csv':
            self.model = pd.read_csv(
                filepath_or_buffer=str(model_path)
                )
        self.max_isos = max(self.model["n_isos"])
        #$fill the data in the rate_header order.  we can reorder if need be
        self.out_path = out_path
        self.graph_folder_isotopes = graph_folder_isotopes
        self.graph_folder_optimization = graph_folder_optimization
        if settings.recognize_available_cores is True:
            self._n_processors = mp.cpu_count()
        else:
            self._n_processors = settings.n_processors
        #$if multiprocessing need to set that up. more than 60 cores causes problems for windows
        if self._n_processors > 60:
            self.n_processors = 60
        self.make_final_header()

    def write(self):
        self.final_rates_model.to_csv(
            path_or_buf=self.out_path,
            index=False
        )
    
    #$this function governs which types of calculations to do and which
    #$ rate equation to use. collects results and puts them in self.rate_model
    def calculate(self):
        start_time = time.time()
        self.initial_filter()
        self.perform_rate_calculations()
        print ("Time Taken: ", time.time() -start_time)
    
    #$make needed header
    def make_final_header(self):
        self.final_header = ["Subject ID", "Protein ID", "Protein Name", 
                             "Sequence", "Abundance rate", "Unique Timepoints",
                             "Number of Measurements",
                             "Approximate Variance",
                             "mean of all absolute residuals",
                             "num times",
                             "n_isos",
                             "num measurements",
                             "time values",
                             "dropped points",
                             "M0 constantly decreasing",
                             "Error column"]
        
        
    #$since the calculation is more complicated than in other modules there are more functions to govern it
    #$this function seperates the data by subject and then calculates the enrichment spline. this is then
    #$passed to run_parallel_groupby which governs the multiprocessing
    def perform_rate_calculations(self):
        #$list to hold output
        rate_calculation_output = []
        #$first we will group by the sample_id so we can grab the water equation 
        #$we need to split by group anyway individual enrichment curves means we need
        #$different calculations for the same peptide in different subjects
        for sample_id, id_df in tqdm(self.model.groupby(group_column1), desc="Subject: "):
            
            max_sample_time = id_df[time_column].max()
            #$need time sorted for the error minimization function.  will sort here to save time
            id_df = id_df.sort_values(by = time_column)
            
            #$since values are the same pulling the top row makes as much sense as any other row
            #$set up the spline
            first_row = id_df.iloc[0].to_dict() 
            subject_x_values = [float(x) for x in first_row["Time Enrichment"].split(", ")]
            subject_y_values = [float(x) for x in first_row["Enrichment Values"].split(", ")]
            spline_object =pickle_safe_spline(subject_x_values, subject_y_values)
                
            #$multiple groupbys should improve readability and is faster than merging them    
            #$will multiprocess the second groupby. since we are grouping by id and sequence
            #$various charge states will be grouped into this data.
            sequence_groupby = id_df.groupby(group_column2)
            one_sample_results =  self.run_parallel_groupby(sequence_groupby, 
                                                           sample_id, 
                                                           spline_object,
                                                           max_sample_time
                                )
            rate_calculation_output.extend(one_sample_results)
        self.final_rates_model = pd.DataFrame(rate_calculation_output, columns = self.final_header)       
                    
    #$run the analysis.  this function doesn't have any calculation itself 
    #$it prepares a function for multiprocessing and thne begins the multiprocessing
    #$uses second answer https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby accessed 12/16/2020
    def run_parallel_groupby(self, groupby_object, sample_id, sample_spline, max_sample_time):
        if settings.debug_level == 0:
            mp_pools = mp.Pool(self._n_processors)
            #$have to pass the settings explicitly (the global settings doesn't play nice with the multiprocessing)
            mp_func = partial(RateCalculator._rate_calculation,
                              sample_id = sample_id,
                              sample_spline =sample_spline,
                              graph_folder_isotopes = self.graph_folder_isotopes,
                              graph_folder_optimization = self.graph_folder_optimization, 
                              minimum_nonzero_points = settings.min_non_zero_timepoints_rate,
                              error_estimation = settings.error_estimation,
                              lowest_allowed_norm_isotope = settings.lowest_allowed_norm_isotope,
                              highest_allowed_norm_isotope = settings.highest_allowed_norm_isotope,
                              m0_decreasing_allowed_noise = settings.m0_decreasing_allowed_noise,
                              median_absolute_residuals_cutoff_single_point = settings.median_absolute_residuals_cutoff_single_point,
                              median_absolute_residuals_cutoff_two_points = settings.median_absolute_residuals_cutoff_two_points,
                              median_absolute_residuals_cutoff_general = settings.median_absolute_residuals_cutoff_general,
                              desired_points_for_optimization_graph = settings.desired_points_for_optimization_graph,
                              max_sample_time = max_sample_time,
                              graph_type = settings.graph_output_format
                              )
            #list is needed to make the tqdm work.  imap_unordered only refers to the order the groupbys entries are analyzed in it maintains the time sort of the data. i did check that.
            results_list = list(tqdm(mp_pools.imap_unordered(mp_func,groupby_object), total = len(groupby_object), desc="Peptide Rate Calculation: "))
            mp_pools.close()
            mp_pools.join()
        else:
            #$no multiprocessing, for troubleshooting
            results_list = []
            mp_func = partial(RateCalculator._rate_calculation,
                              sample_id = sample_id,
                              sample_spline =sample_spline,
                              graph_folder_isotopes = self.graph_folder_isotopes,
                              graph_folder_optimization = self.graph_folder_optimization, 
                              minimum_nonzero_points = settings.min_non_zero_timepoints_rate,
                              error_estimation = settings.error_estimation,
                              lowest_allowed_norm_isotope = settings.lowest_allowed_norm_isotope,
                              highest_allowed_norm_isotope = settings.highest_allowed_norm_isotope,
                              m0_decreasing_allowed_noise = settings.m0_decreasing_allowed_noise,
                              median_absolute_residuals_cutoff_single_point = settings.median_absolute_residuals_cutoff_single_point,
                              median_absolute_residuals_cutoff_two_points = settings.median_absolute_residuals_cutoff_two_points,
                              median_absolute_residuals_cutoff_general = settings.median_absolute_residuals_cutoff_general,
                              desired_points_for_optimization_graph = settings.desired_points_for_optimization_graph,
                              max_sample_time = max_sample_time,
                              graph_type = settings.graph_output_format
                              )
            for g in groupby_object:
                results_list.append(mp_func(g))
        
        return results_list
    
    
    #$for now this drops this just drops error messages in certain columns
    def initial_filter(self):
        self.model["Theoretical Unlabeled Normalized Abundances"] = \
            self.model["Theoretical Unlabeled Normalized Abundances"].apply(
                RateCalculator._check_for_error_message)
        
        self.model = self.model.dropna(subset = ["Theoretical Unlabeled Normalized Abundances"])    
    
    @staticmethod
    def _check_for_error_message(value):
        try:
            test = value.split(", ")
            test = [float(x) for x in test]
            return test
        except:
            return np.nan
    
    #$this is the return function for _rate_calculation
    #$since it is a pandas series and there are many values, this function was made to make the 
    #$_rate_calculation function more readable since it can error out and need output at several different locations
    @staticmethod
    def _make_return_series(sample_id, protein_name, protein_english_name,
                                 sequence_name, rate, time_data, unique_times,n_isos,
                                 dropped_points, error, m0_constantly_decreasing, 
                                 residual_error, approx_variance):
        return_series = pd.Series({
                "Subject ID": sample_id,
                "Protein ID": protein_name,
                "Protein Name": protein_english_name,
                "Sequence": sequence_name,
                "Abundance rate": rate,
                "Unique Timepoints": len(unique_times), 
                "Number of Measurements": len(time_data),
                "Approximate Variance": approx_variance,
                "mean of all absolute residuals":residual_error,
                "num times": len(time_data),
                "n_isos": n_isos,
                "num measurements": len(time_data) *n_isos,
                "time values": "; ".join(str(x) for x in time_data),
                "dropped points": dropped_points,
                "M0 constantly decreasing": m0_constantly_decreasing,
                "Error column":error
                })
        return return_series

    #$perform the rate calculation
    @staticmethod
    def _rate_calculation(groupby_tuple, sample_id, sample_spline, 
                          graph_folder_isotopes, graph_folder_optimization,
                          minimum_nonzero_points, error_estimation,
                          lowest_allowed_norm_isotope,
                          highest_allowed_norm_isotope,
                          m0_decreasing_allowed_noise,
                          median_absolute_residuals_cutoff_single_point,
                          median_absolute_residuals_cutoff_two_points,
                          median_absolute_residuals_cutoff_general,
                          desired_points_for_optimization_graph, 
                          max_sample_time,
                          graph_type):

        sequence_name = groupby_tuple[0]
        seq_df = groupby_tuple[1]
        
        
        
        #$ pull out the first row to grab things that are consistent, like n value and n isos
        generic_first_row = seq_df.iloc[0].to_dict() #same logic as above
        n_isos = int(generic_first_row["n_isos"])
        n_value = round(generic_first_row["literature_n"])
        protein_name = generic_first_row["Protein ID"]
        protein_english_name =  generic_first_row["Protein Name"]
        time_data = np.asarray(seq_df[time_column])
        #$get the data out 
        empirical_data = np.asarray(seq_df["abundances"])
        
        #$confirm the data is good before complex calculations begin
        #$ currently confirms that there are as many non-zero timepoints as the user specified.
        unique_times = sorted(list(set(time_data)))
        non_zero_times = [t for t in unique_times if t!=0.0]
        if len(non_zero_times) < minimum_nonzero_points:
            return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, "", time_data, unique_times,n_isos, "",
                                 "insufficient non-zero timepoints", "", "", "")
            return return_series
        
        #$theoretical y0 from emass
        y0 = generic_first_row["Theoretical Unlabeled Normalized Abundances"]
        
        #$need to collect the data from their strings and perform the normalization
        normed_isotope_data= []
        for measured_isotopes in empirical_data:
            measured_isotopes = measured_isotopes[1:-1] #$trim terminal parentheses
            #$make an array of floats from the measured_isotopes string
            measured_isotopes = np.asarray([float(x) for x in measured_isotopes.split(", ")])
            #$perform the normalization
            normed_isotope_data.append(normalize(measured_isotopes))
            
        normed_isotope_data = np.asarray(normed_isotope_data[:len(time_data)])
        full_length = len(normed_isotope_data)
        
        #$we need to trim out points that are theoretically impossible.
        #$we'll do that by calculating m0 peaks going up or m1-4 peaks going down by more than a certain amount
        #$the certain amounts will be based on their abundance with no extra isotope
        #$allowed error can be adjusted in the settings
        y0_allowed_error = [y0[0] * highest_allowed_norm_isotope]
        for i in range(1,n_isos):
            y0_allowed_error.append(y0[i] * lowest_allowed_norm_isotope)
        row_mask = (normed_isotope_data[:, 0] < y0_allowed_error[0]) & (normed_isotope_data[:, 1:] > y0_allowed_error[1:]).all(1) 
        
        normed_isotope_data = normed_isotope_data[row_mask]
        time_data = time_data[row_mask]
        #$recalculate unique_times now that we may have dropped some of them
        unique_times = sorted(list(set(time_data)))
        
        non_zero_times = [t for t in unique_times if t!=0.0]
        dropped_points = full_length - len(normed_isotope_data)
        
        
        #$since we just fitered out theoretically impossible points we need to ensure we still
        #$have sufficient non-zero points to perform the fit
        if len(non_zero_times) < minimum_nonzero_points:
            return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, "", time_data, unique_times,n_isos,
                                 dropped_points,
                                 "insufficient non-zero timepoints after filtering",
                                 "", "", "")
            return return_series
        
        #$warn if m0 is notstrictly decreasing only functions as a filter if there is only one point per time
        #$this is to catch fits with a large dip or spike in data that results in a terrible fit.  a noise increase
        #$should not trigger this. usually only affects one point so multiple points (charge states) at each time are
        #$resistant to this so that reports an error but will not drop
        
        #$unfortunately the possibility for multiple timepoints complicates the analysis since order the points are analyzed
        m0_values = normed_isotope_data[:, 0]
        relative_m0_decreasing_allowed_noise = y0[0] *m0_decreasing_allowed_noise
        if len(unique_times) == len(time_data):
            for r in range(1, len(m0_values)):
                if m0_values[r] > (m0_values[r-1] + relative_m0_decreasing_allowed_noise):
                    return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, "", time_data, unique_times,n_isos,
                                 dropped_points,
                                 "m0 is not constantly decreasing when there is only one point per time",
                                 "", "", "")
                    return return_series
            else:
                m0_constantly_decreasing = True
        else:
            needed_dict = defaultdict(list)
            for k,v in zip(time_data, m0_values):
                needed_dict[k].append(v)
            for t in range(1, len(unique_times)):
                if max(needed_dict[unique_times[t]]) > (min(needed_dict[unique_times[t-1]]) + relative_m0_decreasing_allowed_noise):
                    m0_constantly_decreasing = False
                    break
            else:
                m0_constantly_decreasing = True
        
        #$these functions used for the fit must be defined here.  They were developed by Christian in the Transtrum lab
       
        def dydt(t,y,k,ibase,nobs): #vector of velocity towards equilibrium
            rho=[binom.pmf(ix,n_value,sample_spline.final_spline(t)) for ix in range(nobs)] #equilibrium distribution for the current deut-level
            rho=np.convolve(ibase,rho)[0:nobs]
            rho=rho/sum(rho)
            return [rhoi*(k)-k*yi for (rhoi,yi) in zip(rho,y)] #the movement towards that equilibrium
        
        def mlevels2(t,k,ibase,nobs): #return simulated M-distributions at all times t for a given k
            mval=scipy.integrate.odeint(dydt,ibase+sample_spline.final_spline(0),t,args=tuple([k,ibase,nobs]),tfirst=True)
            #$quick normalization or they will all be slightly off
            sum_array = mval.sum(axis = 1)
            return np.divide(mval.T,sum_array).T
        
        def sse(k,t,y,ibase,nobs):#sum squared error
            return sum((mlevels2(t,k,ibase,nobs) - y).flatten()**2)
        
        number_of_non_zero_fit_points = len([t for t in time_data if t != 0 ])
        #$the 0 in 0:n_isos used to be 1.  unsure why as it was cutting data.  this may have been for removal of 0 in testing
        #$will keep the 0:n_isos form as a legacy of this to aid in potential troubleshooting
        try:
            #$if we don't allow 0 non-zero timepoints we don't need to check here and can assume that there is at least 1 non-zero timepoint 
            #$the fitting equation assumes that a 0 is present.  however 0 is forced by the theoretical value
            #$(like how y= mx has to go through 0).  therefore if 0 is not there we'll add it.  the perfect 0 won't skew the fit, it will prevent
            #$the first data point from effectively being thrown out
            if 0 not in unique_times:
                fitting_times = np.insert(time_data, 0, 0)
                fitting_values = [np.asarray(y0)]
                for l in range(len(normed_isotope_data)):
                    fitting_values.append(normed_isotope_data[l])
                fitting_values = np.asarray(fitting_values)
                added_zero = True #$need to know for residual checks
            else:
                fitting_times = time_data
                fitting_values = normed_isotope_data
                added_zero = False
            #$perform the error minimization fit with the ODE functions
            bmin=brent(sse,args=(fitting_times,fitting_values,y0, n_isos),brack=(1e-6,1e-4))
            k_value = bmin
        #$ if there was a fit warning treat it as an error
        except ODEintWarning:
            return_series = RateCalculator._make_return_series(sample_id, 
                             protein_name, protein_english_name,
                             sequence_name, "", time_data, 
                             unique_times,n_isos,dropped_points,
                             "Fit is too poor for accurate calculation",
                             m0_constantly_decreasing, "", "")
            return return_series
        
        #$we will use the mean of the absolute residuals as an error metric so we need to calculate that
        #$first we need to determine what the fit predicts the points are
        calculated_values_at_each_time = mlevels2(fitting_times, k=k_value, ibase = y0, nobs=n_isos)
        #$same problem as above. can easily use fitting times to avoid, but then need to cut out any zero we addedwhen calculating residuals
        if added_zero:
            delta_rel_isotope_int = normed_isotope_data- calculated_values_at_each_time[1:]
        else:
            delta_rel_isotope_int = normed_isotope_data- calculated_values_at_each_time
        #$absolute value the error to prevent the error in opposite directions from cancelling out and giving an overly optomistic view
        abs_delta = abs(delta_rel_isotope_int)
        
        #$metric 1 average of the abs value of the errors
        absolute_error_total_mean = np.mean(abs_delta)
        
        #$approximate fit variation. used as an error metric and is therefore outside of the graphers
        try:
            h=bmin/7.38906 #h is 2 log-units less than the min
            xval=bmin+[-h,0,h]
            yfit=[0,0,0]
            #$if we try single point with just one point here it sends the error massive positive or negative
            #$and even on multiple points sse needs a zero value so this should help
            for j,x in enumerate(xval):
                yfit[j]=sse(x,fitting_times,fitting_values,y0, n_isos)
            parfit,pcov=curve_fit(parabola, xval, yfit)
            approx_variance = .5/parfit[0] #f''=2a for parabola, var = 1/f''
        except ODEintWarning:
            approx_variance = "could not find the variance"
        
        
        #$now that we have the data and fit we need we can start making graphs
        
        #$first we need the theoretically perfect lines so we can grph them
        theory_times = np.arange(0, max_sample_time + 1, 0.5)
        predicted_isotope_values = mlevels2(theory_times, k=k_value, ibase=y0, nobs = n_isos)
        
        graph_file_name = sample_id + "_" + sequence_name + "_isotopes."+ graph_type
        graph_name_isotopes = os.path.join(graph_folder_isotopes, graph_file_name )
        if error_estimation != "none":
            graph_file_name_optimize = sample_id + "_" + sequence_name + "_optimization." +graph_type
            graph_name_optimize = os.path.join(graph_folder_optimization, graph_file_name_optimize)
        
        #$now let's actually graphs. rate first
        graph_rate_results(n_isos, graph_name_isotopes, time_data, 
                       normed_isotope_data, predicted_isotope_values, theory_times,
                       y0,
                       sample_id + "_" + sequence_name, graph_type)
        #$the user specifies the way they want the error represented or even if they want it presented at all.
        #$this if takes care of that.
        if error_estimation == "approximate":
            error_xv=np.arange(0,bmin*3,bmin*3/desired_points_for_optimization_graph)
            error_array = parabola(error_xv,parfit[0],parfit[1],parfit[2])
            legend_name = "approx for f''"
            graph_optimization_of_error(k_value, error_xv, error_array, graph_name_optimize, sample_id + "_" + sequence_name, legend_name, graph_type)
        elif error_estimation == "exact":
            error_xv=np.arange(0,bmin*3,bmin*3/desired_points_for_optimization_graph)
            legend_name = "real cost of k"
            #$possible to error out here so catch it if so
            try :
                #$mimic up top.
                if 0 not in time_data:
                    time_data = np.insert(time_data, 0,0)
                    normed_isotope_data = np.insert(normed_isotope_data, 0, y0,0)
                error_array=[sse(x,time_data,normed_isotope_data,y0, n_isos) for x in error_xv]
                graph_optimization_of_error(k_value, error_xv, error_array, graph_name_optimize, sample_id + "_" + sequence_name, legend_name, graph_type)
            except ODEintWarning:
                error_array = "Error"
        
        #$using single or two point fits generally need tighter filters, so for this evaluation we will employ those filters here
        if number_of_non_zero_fit_points == 1 and absolute_error_total_mean > median_absolute_residuals_cutoff_single_point:
            return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, k_value, time_data, 
                                 unique_times,n_isos, dropped_points,
                                 "mean of the absolute residuals is too high", 
                                 m0_constantly_decreasing, absolute_error_total_mean, approx_variance)
            return return_series
        elif number_of_non_zero_fit_points == 2 and absolute_error_total_mean > median_absolute_residuals_cutoff_two_points:
            return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, k_value, time_data, 
                                 unique_times,n_isos, dropped_points,
                                 "mean of the absolute residuals is too high", 
                                 m0_constantly_decreasing, absolute_error_total_mean, approx_variance)
            return return_series
        elif number_of_non_zero_fit_points > 2 and absolute_error_total_mean > median_absolute_residuals_cutoff_general:
            return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, k_value, time_data, 
                                 unique_times,n_isos, dropped_points,
                                 "mean of the absolute residuals is too high", 
                                 m0_constantly_decreasing, absolute_error_total_mean, approx_variance)
            return return_series
        
        #$if all is good this is the final row calculated
        return_series = RateCalculator._make_return_series(sample_id, 
                                 protein_name, protein_english_name,
                                 sequence_name, k_value, time_data, 
                                 unique_times,n_isos, dropped_points,
                                 "", m0_constantly_decreasing, 
                                 absolute_error_total_mean, approx_variance)
        return return_series