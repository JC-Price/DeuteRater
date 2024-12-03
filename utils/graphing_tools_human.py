# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2021 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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

"""
functions for making the various graphs
"""

error_line_symbol = 'k--'
data_points_symbol = 'ro'

multi_colors = ['r', 'g', 'b', 'c', 'm', 'y']


# need to graph the results of the binomial fit
# adapted from code by Christian in the Transtrum lab
def graph_rate_results(n_isos, save_file_name, time_data,
                       normed_isotope_data, mval, theory_times, theory_zero,
                       subject_sequence, save_format):
    plt.title(f"Isotope levels {subject_sequence}")
    plt.xlabel("Time")
    plt.ylabel("Fraction")
    for j in range(n_isos):
        plt.plot(theory_times, mval[:, j], multi_colors[j] + '-', label="M" + str(j))
    for j in range(n_isos):
        plt.plot(time_data, normed_isotope_data[:, j], multi_colors[j] + "o")
        plt.plot(0, theory_zero[j], multi_colors[j] + "d")
    plt.legend(loc="upper right")
    plt.savefig(save_file_name, format=save_format)
    plt.clf()


# need to graph the optimization of error for troubleshooting purposes
# adapted from code by Christian in the Transtrum lab
def graph_optimization_of_error(k_value, theoretical_k, cost, save_file_name, subject_sequence, legend_name,
                                save_format):
    plt.title(f"Finding optimal k {subject_sequence}")
    plt.ylabel("Sum-square Error")
    plt.xlabel("k")
    plt.plot(theoretical_k, cost, "bo-", label=legend_name)
    plt.plot(np.repeat(k_value, 2), [min(cost), max(cost)])
    plt.legend()
    plt.savefig(save_file_name, format=save_format)
    plt.clf()


# graph the enrichment spline
def enrichment_graph(x_values, y_values, predicted_x,
                     predicted_y, subject_name, save_file, save_format):
    plt.title(subject_name)
    plt.xlabel("Time")
    plt.ylabel("Enrichment")
    plt.plot(x_values, y_values, data_points_symbol)
    plt.plot(predicted_x, predicted_y, error_line_symbol)

    plt.savefig(save_file, format=save_format)
    plt.clf()


# graph the protein roll up
def graph_average_boxplots(values, save_file_name, subject_protein_name, save_format):
    # give the values some jitter, so they are easily visible from each other
    # jitter is entirely random.  if we decide to remove, just make all  x_values 1
    x_values = np.random.normal(1, .035, size=len(values))
    average = np.mean(values)

    plt.plot(x_values, values, 'r*', markersize=11, label="Peptide Rates")
    # since the line width is an argument, is needs a value, but if we try and declare it as a dict,
    # it will not be defined.  declaring its value and forcing to dict gets around that
    plt.boxplot(values,
                whiskerprops=dict(linewidth=1.5),
                boxprops=dict(linewidth=3.0),
                capprops=dict(linewidth=2.0),
                medianprops=dict(linewidth=2.0)
                )
    plt.plot([1.00], [average],
             markerfacecolor="deepskyblue",
             markeredgecolor="deepskyblue",
             markersize=10,
             marker="d",
             label="Mean",
             linestyle="None")
    plt.legend()
    plt.title(subject_protein_name)
    plt.xticks([])  # remove x-axis labels to avoid confusion
    plt.ylabel('Sequence Rate', fontsize=12)
    plt.savefig(save_file_name, format=save_format)
    plt.clf()
