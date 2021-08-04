# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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
import os

main_line_symbol = 'k-'
#error_line_symbol = 'k--'
data_points_symbol = 'ro'
MAXIMUM_GRAPH_RATE_ERROR = 5
MINIMUM_GRAPH_RATE_ERROR = -2

#$the colon is not allowed but seems to make an empty file with a partial name.
#$either way the check is here to prevent problems if it is necessary
bad_save_file_characters = ["/", "\\", ":", "*", "?", "\"", "<", ">", "|"]

def graph_rate(name, x_values, y_values, rate, asymptote, ci, rate_equation, 
               save_folder_name, maximum, asymptote_option, biomolecule, errors = [], full_data=None, title=None):
    #$start naming things
    if title is None:
        title = name
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel(f'Fraction of Total {biomolecule} that is new')
    
    #$make the values for the line
    #$make_error_lines is just to remove unoptimal fits for plus and minus
    #$the main should have triggered the try except in rate_calculator.py if
    #$it is the issue.
    make_error_lines = True
    fit_line_x = np.arange(0,maximum + maximum/10, .1)
    if asymptote_option == 'variable' or asymptote_option == "Variable":
        fit_line_y = rate_equation(fit_line_x, k = rate, a  = asymptote)
        try:
            fit_line_y_minus_error = rate_equation(fit_line_x, 
                                               k = rate - ci, a  = asymptote)
            fit_line_y_plus_error =  rate_equation(fit_line_x, 
                                               k = rate + ci, a  = asymptote)
        except:
            make_error_lines = False
    else:
        fit_line_y = rate_equation(fit_line_x, rate)
        try:
            fit_line_y_minus_error = rate_equation(fit_line_x, rate - ci)
            fit_line_y_plus_error = rate_equation(fit_line_x, rate + ci)
        except:
            make_error_lines = False
        
    #$ if add multiple conditions put each in parentheses  and use | or &
    #$ unfortunately nans can cause errors.  only seem to occur if using one timepoint and only error in an .exe but may as well sort this out
    plt.plot(fit_line_x, fit_line_y, main_line_symbol)
    if make_error_lines and not np.isnan(fit_line_y_plus_error).any() and not np.isnan(fit_line_y_minus_error).any():
        fit_line_y_plus_error[fit_line_y_plus_error > MAXIMUM_GRAPH_RATE_ERROR] = \
            MAXIMUM_GRAPH_RATE_ERROR
        fit_line_y_minus_error[fit_line_y_minus_error < MINIMUM_GRAPH_RATE_ERROR] = \
            MINIMUM_GRAPH_RATE_ERROR   
    
    #$plot  lines and points
    if full_data is None:
        plt.plot(x_values, y_values, data_points_symbol)
    else:
        color_list = {"M+H": 'r',
                      "M+Na": 'b',
                      "M+NH4": 'g',
                      'M+H-[H20]': 'magenta',
                      "M+Na-[H2O]": 'gold',
                      "M+NH4-[H2O]": 'lime'}
        shapes = {"M1": 'o', "M2": 's', "M4": '^'}
        charge_fill = {1: False, 2: True}
        adduct_groups = full_data.groupby("Adduct")
        for adduct_group in adduct_groups:
            adduct = adduct_group[0]
            color = color_list[adduct]
            adduct_data = adduct_group[1]
            rep_groups = adduct_data.groupby("bio_rep")
            for rep_group in rep_groups:
                rep = rep_group[0]
                shape = shapes[rep]
                rep_data = rep_group[1]
                charge_groups = rep_data.groupby("z")
                for charge_group in charge_groups:
                    charge = charge_group[0]
                    should_fill = charge_fill[charge]
                    charge_data = charge_group[1]
                    n_val = charge_data['n_value'].unique()[0]
                    x = charge_data['time']
                    y = charge_data['abund_fn']
                    label = f'{rep}_z{charge}_{adduct} n={n_val}'
                    if should_fill:
                        plt.scatter(x, y, marker=shape,
                                    facecolor=color, edgecolor='k',
                                    label=label)
                    else:
                        plt.scatter(x, y, marker=shape,
                                    facecolor='none', edgecolor=color,
                                    label=label)
    if make_error_lines:
        plt.fill_between(fit_line_x, fit_line_y_minus_error, fit_line_y_plus_error, color = 'black', alpha = .15)
        #plt.plot(fit_line_x, fit_line_y_plus_error, error_line_symbol)
        #plt.plot(fit_line_x, fit_line_y_minus_error, error_line_symbol)
    if errors != [] : #$ only if roll up so need error bars
        plt.errorbar(x_values, y_values, yerr = errors,  elinewidth = 1, 
            ecolor = 'red', linewidth = 0)
    #$save figure and clear it for next time
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_size("small")
    plt.legend(prop=font)
    try:
        filename = os.path.join(save_folder_name, name) + ".png"
        
        
        plt.savefig(filename)

    #$the following characters are not allowed in windows file names: /\ : * ? " <> |
    #$replace them with underscores (yes this could erase duplicates, but this is unlikely 
    #$first place, that situation is unlikely problematic enough to bother with)
    except OSError:
        for bad_char in bad_save_file_characters:
           name = name.replace(bad_char, "_")
        filename = os.path.join(save_folder_name, name) + ".png"
        plt.savefig(filename)

    plt.clf()

