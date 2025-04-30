# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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

deal with splines and fitting the equations
will also do any error checking if necessary 

"""

import os
import numpy as np
import pandas as pd
import scipy.interpolate as si

from utils.graphing_tools_human import enrichment_graph

# the spline equation will be separate for easier calling in other modules
def spline_interp(t, d):
    std_d = np.std(d)
    w_y = [3/std_d for i in range(len(t))]
    spl = si.splrep(t, d, w=w_y)
    return lambda x: si.splev(x, spl)

# this just holds some data for ease
class subject_data(object):
    def __init__(self, subject_name, x, y):
        self.subject_name = subject_name
        self.x = x
        self.y = y
    def calculate_spline(self):
        self.spline =  spline_interp(self.x, self.y)
    
    def calculate_theory_data(self):
        self.theory_x = np.arange(0, max(self.x), .1)
        self.theory_y = self.spline(self.theory_x)

old_header =["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
# we'll use a class for consistency with other parts
class PerformEnrichmentClass(object):
    def __init__(self, input_filename, graph_folder_name, graph_type):
        temp_df = pd.read_csv(input_filename, sep = "\t")
        self.subject_ids = list(temp_df[old_header[0]].unique())
        # if the number of enrichment measurments is less than the number of subject samples we can end with nans in the subject list
        # we'll remove that.
        self.subject_ids = [s for s in self.subject_ids if pd.notnull(s)]
        self.load_previous_data(input_filename)
        self.graph_folder_name = graph_folder_name
        self.graph_type = graph_type
        
    def perform_calculations(self):
        # sometimes a fit fails.  instead of raising a warning or error, the spline just returns nans
        # fortunately the spline can extrapolate and should start at or near 0.  so to prevent errors in rate calculation, we'll try it,
        # track it, and then return it later
        self.error_list = []
        for subject_key in self.subject_dictionary.keys():
            temp_subject_class = self.subject_dictionary[subject_key]
            temp_subject_class.calculate_spline()
            temp_subject_class.calculate_theory_data()
            temp_graph_name = os.path.join(self.graph_folder_name, temp_subject_class.subject_name +f".{self.graph_type}")
            enrichment_graph(temp_subject_class.x, temp_subject_class.y, 
                             temp_subject_class.theory_x, temp_subject_class.theory_y,
                             temp_subject_class.subject_name, temp_graph_name, self.graph_type)
            #  report error might as well graph it anyway.  
            if any(np.isnan(temp_subject_class.theory_y)):
                self.error_list.append(subject_key)

    def report_error(self):
        if self.error_list == []:
            return ""
        error_message = "The following Subjects failed to fit: \n{}\nlikely because of errors in entering data in the enrichment table. Analysis has stopped. Please try again.".format("\n".join(self.error_list))
        return error_message
    
    # need to load the file and prepare for fitting
    def load_previous_data(self, input_filename):
        self.subject_dictionary = {}
        df = pd.read_csv(input_filename, delimiter=("\t"))
        df = df[old_header]
        for subject in self.subject_ids:   
            temp_df = df[df[old_header[0]] == subject]
            #  the spline needs to be in order. this is faster and easier than trying to do it when making the output
            temp_df = temp_df.sort_values(old_header[1])
            # need to pass numpy array values for the fits
            self.subject_dictionary[subject] = subject_data(
                subject,
                temp_df[old_header[1]].to_numpy(),
                temp_df[old_header[2]].to_numpy()
                )


if __name__ == '__main__':
    test = PerformEnrichmentClass("", "")
    