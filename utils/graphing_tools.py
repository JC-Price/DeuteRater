# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:10:22 2020

@author: JCPrice

holds  basic functions for graphing. should be called as needed, no need to 
run as the main.
we'll call directly from the calculator, which will simplify the code here
by a significant amount
"""

import matplotlib.pyplot as plt
import numpy as np
import os

main_line_symbol = 'k-'
error_line_symbol = 'k--'
data_points_symbol = 'ro'
MAXIMUM_GRAPH_RATE_ERROR = 10
MINIMUM_GRAPH_RATE_ERROR = -5


def graph_rate(name, x_values, y_values, rate, asymptote, ci, rate_equation, 
               save_folder_name, maximum, asymptote_option, errors = []):
    #$start naming things
    plt.title(name)
    plt.xlabel('Time')
    plt.ylabel('Fraction of Total Protein that is new')
    
    
    #$make the values for the line
    #$make_error_lines is just to remove unoptimal fits for plus and minus
    #$the main should have triggered the try except in rate_calculator.py if
    #$it is the issue.
    make_error_lines = True
    fit_line_x = np.arange(0,maximum + maximum/10, .1)
    if asymptote_option == 'variable':
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
    if make_error_lines:
        fit_line_y_plus_error[fit_line_y_plus_error > MAXIMUM_GRAPH_RATE_ERROR] = \
            MAXIMUM_GRAPH_RATE_ERROR
        fit_line_y_minus_error[fit_line_y_minus_error < MINIMUM_GRAPH_RATE_ERROR] = \
            MINIMUM_GRAPH_RATE_ERROR   
    
    #$plot  lines and points
    plt.plot(fit_line_x, fit_line_y, main_line_symbol)
    if make_error_lines:
        plt.plot(fit_line_x, fit_line_y_plus_error, error_line_symbol)
        plt.plot(fit_line_x, fit_line_y_minus_error, error_line_symbol)
    plt.plot(x_values, y_values, data_points_symbol)
    if errors != [] : #$ only if roll up so need error bars
        plt.errorbar(x_values, y_values, yerr = errors,  elinewidth = 1, 
            ecolor = 'red', linewidth = 0)
    #$save figure and clear it for next time
    filename = os.path.join(save_folder_name, name) + ".pdf"
    plt.savefig(filename)
    plt.clf()

