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
Need to merge the sequences into proteins.
most of this will consist of various filters remove problematic rates
will set up to allow multiprocessing. it is likley that this will be unnecessary
initially because it is just filtering, but it will allow flexibility later
also this will not be the rate limiting step so multiprocessing overhead is of
limited concern

for right now we will keep all of the filters here.  there is an argument to be 
made for doing the filtering earlier, as the graphs are made and just merging
this with the previous step. However, that needs to complicate the moste complicated
class, and add sevral booleans for things like if we should make graphs or not.
Therefore we will start here, and adjust later as needed.
"""
import pandas as pd
import multiprocessing as mp
import numpy as np
import os

from functools import partial
from pathlib import Path
from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
from utils.graphing_tools_human import graph_average_boxplots

# will be grouping by subject and protein id
group_column1 = "Subject ID"
group_column2 = "Protein ID"
rate_column = "Abundance rate"
error_column = "Error column"
# all text errors that rate_calculator.py can report for a fit.
# if any new errors are created add them here.
text_errors = ["insufficient non-zero timepoints", 
               "insufficient non-zero timepoints after filtering", 
               "Fit is too poor for accurate calculation",
               "m0 is not constantly decreasing when there is only one point per time",
               "mean of the absolute residuals is too high"]

# as with all the calculation steps this is a class for consistent calls in the main
class Peptides_to_Proteins(object):
    def __init__(self, model_path, out_path, settings_path, graph_folder_path):
        settings.load(settings_path)
        model_path = Path(model_path)
        # allow .tsv or .csv files
        if model_path.suffix == '.tsv':
            self.model = pd.read_csv(
                filepath_or_buffer=str(model_path),
                sep='\t'
                )
        elif model_path.suffix == '.csv':
            self.model = pd.read_csv(
                filepath_or_buffer=str(model_path)
                )
        self.out_path = out_path
        # if multiprocessing need to set that up. more than 60 cores causes problems for windows
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.80)
        else:
            self._n_processors = settings.n_processors
        if self._n_processors > 60:
            self.n_processors = 60
        
        self.graph_folder = graph_folder_path
        self.make_final_header()
    
    def write(self):
        self.final_rates_model.to_csv(
            path_or_buf=self.out_path,
            index=False
        )
    
    # make the final header and add to the error headers
    def make_final_header(self):
        self.final_header = ["Subject ID", "Protein ID", "Protein Name", "Num Sequences Used",
                             f"{settings.protein_combination_method} Rate", "Rate Std Dev", 
                             "Num Peptides Above Max. Allowed Rate", 
                             "Num Peptides Below Min. Allowed Rate"]
        for text_error in text_errors:
            self.final_header.append("number of sequences with error \"{}\"".format(text_error))
     # run the analysis.  this function doesn't have any calculation itself (other than merging results and header adjustment)
    # it prepares a function for multiprocessing and thne begins the multiprocessing
    def calculate(self):
        self.make_final_header()
        if settings.debug_level == 0:
            mp_pools = mp.Pool(self._n_processors)
            
            # unlike rate_calculator.py we'll groupby both at once.  if we need to do anything within the 
            # sample groups split it out like rate calculator
            groupby_object = self.model.groupby([group_column1, group_column2])
            mp_func = partial(Peptides_to_Proteins._mp_grouping_function,
                              minimum_rate = settings.minimum_allowed_sequence_rate,
                              max_rate = settings.maximum_allowed_sequence_rate,
                              minimum_required_rates = settings.minimum_sequences_to_combine_for_protein_rate,
                              text_errors = text_errors,
                              graph_folder = self.graph_folder,
                              graph_type = settings.graph_output_format,
                              roll_up_type = settings.protein_combination_method)
                              
            results_list = list(tqdm(mp_pools.map(mp_func,[group_tuple for group_tuple in groupby_object]), 
                                     total = len(groupby_object), desc = "Protein Roll Up: "))
            mp_pools.close()
            mp_pools.join()
        
        else:
            # no multiprocessing for troubleshooting 
            results_list = []
            groupby_object = self.model.groupby([group_column1, group_column2])
            mp_func = partial(Peptides_to_Proteins._mp_grouping_function,
                              minimum_rate = settings.minimum_allowed_sequence_rate,
                              max_rate = settings.maximum_allowed_sequence_rate,
                              minimum_required_rates = settings.minimum_sequences_to_combine_for_protein_rate,
                              text_errors = text_errors)
            for g in groupby_object:
                results_list.append(mp_func(g))
      
        self.final_rates_model = pd.DataFrame(results_list, columns = self.final_header)   

        # simplify at the end.  makes troubleshooting and adjustment easier
        # and should not drastically increase the time taken
        if not settings.verbose_output:
            self.trim_verbose_data()        
        
    # if the user doesn't want the number and types of errors they can be left out for ease of viewing
    def trim_verbose_data(self):
        needed_columns = ["Subject ID", "Protein ID", "Protein Name", "Num Sequences Used", f"{settings.protein_combination_method} Rate", "Rate Std Dev"]
        self.final_rates_model = self.final_rates_model[needed_columns] 
       
    #  the actual function to do the caluclation
    @staticmethod
    def _mp_grouping_function(groupby_element, minimum_rate, max_rate,
                              minimum_required_rates, text_errors, graph_folder,
                              graph_type, roll_up_type):
        id_tuple, id_dataframe = groupby_element
        return_dict = {"Subject ID": id_tuple[0], "Protein ID": id_tuple[1]}
        return_dict["Protein Name"] = list(id_dataframe["Protein Name"])[0]
        
        
        # we're going to apply the filters one at a time so we can collect numbers
        # filtered out. can merge later
        
        # first we need to get rid any non numerical values so we can do math
        for text_error in text_errors:
            start_len = len(id_dataframe.index)
            id_dataframe = id_dataframe[id_dataframe[error_column] != text_error]
            return_dict["number of sequences with error \"{}\"".format(text_error)] = start_len-len(id_dataframe.index)
         
        # the later math checks can run into issues if the data frame is empty (particularly indexing)
        # may shift out of this later if necessary
        if id_dataframe.empty:
            return_dict["Num Peptides Above Max. Allowed Rate"] = 0
            return_dict["Num Peptides Below Min. Allowed Rate"] = 0
            return_dict["Num Sequences Used"] = 0
            return_dict[f"{roll_up_type} Rate"] = "No Rates were calculated for this ID"
            return_dict["Rate Std Dev"] = "No Rates were calculated for this ID"
            return pd.Series(return_dict)
        
        # now that we have dropped all known text errors we can force a numerical conversion
        id_dataframe[rate_column] = pd.to_numeric(id_dataframe[rate_column])
        # will not remove nans since all non-numerical data has been removed.  
        # if there is text still there we need to know it
        
        
        # collect the errors
        num_numerical_values = len(id_dataframe.index)
        id_dataframe = id_dataframe[id_dataframe[rate_column] < max_rate]
        num_not_too_big_values = len(id_dataframe.index)
        id_dataframe = id_dataframe[id_dataframe[rate_column] > minimum_rate]
        num_values_in_allowed_range = len(id_dataframe.index)
        
        return_dict["Num Peptides Above Max. Allowed Rate"] = num_numerical_values - num_not_too_big_values
        return_dict["Num Peptides Below Min. Allowed Rate"] = num_not_too_big_values - num_values_in_allowed_range
        return_dict["Num Sequences Used"] = num_values_in_allowed_range
        
        # determine if we have enough data for calculation we can do our median or average
        if num_values_in_allowed_range == 0:
            return_dict[f"{roll_up_type} Rate"] = "No Good Rates were calculated for this ID"
            return_dict["Rate Std Dev"] = "No Good Rates were calculated for this ID"
        elif num_values_in_allowed_range < minimum_required_rates:
            return_dict[f"{roll_up_type} Rate"] = "insufficient rates to combine"
            return_dict["Rate Std Dev"] = "insufficient rates to combine"
        else:
            if roll_up_type == "Average":
                return_dict[f"{roll_up_type} Rate"] = np.mean(id_dataframe[rate_column])
            elif roll_up_type == "Median":
                return_dict[f"{roll_up_type} Rate"] = np.median(id_dataframe[rate_column])
            return_dict["Rate Std Dev"] = np.std(id_dataframe[rate_column])
        
            protein_graph_title = f"{id_tuple[0]}_{id_tuple[1]}"
            full_graph_name = os.path.join(graph_folder, f"{protein_graph_title}.{graph_type}")
            graph_average_boxplots(id_dataframe[rate_column], full_graph_name, protein_graph_title, graph_type)
            
        return pd.Series(return_dict)