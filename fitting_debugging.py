from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import t

import utils.rate_equations as dur

MAXIMUM_GRAPH_RATE_ERROR = 5
MINIMUM_GRAPH_RATE_ERROR = -2

# x_file = "D:\\DR Testing\\Graphing\\x_coords.csv"
# y_file = "D:\\DR Testing\\Graphing\\y_coords.csv"
x_file = "D:\\DR Testing\\Abundance\\x_coords_ab.csv"
y_file = "D:\\DR Testing\\Abundance\\y_coords_ab.csv"
p0 = [.1, 1]


def main():
    x_coords = pd.read_csv(x_file)
    y_coords = pd.read_csv(y_file)
    
    x = np.concatenate(x_coords.to_numpy())
    y = np.concatenate(y_coords.to_numpy())
    
    rate_eq = partial(dur.simple, p_adj=0.0)
    popt, pcov = curve_fit(f=rate_eq, xdata=x, ydata=y, p0=p0)
    rate = popt[0]
    asymptote = popt[1] if len(popt) > 1 else 1
    confint = t.ppf(.975, 6 - len(popt)) * np.sqrt(np.diag(pcov))[0]
    y_predicted = dur.simple(x, rate, asymptote, 0.0)
    r_2 = dur.calculate_r2(y, y_predicted)
    
    plt.xlabel('Time (Days)')
    plt.ylabel(f'Fraction New')
    
    for time in x:
        plt.axvline(x=time, linewidth=.5, color="k")
    
    maximum = 32.0
    fit_line_x = np.arange(0, maximum + maximum / 10, .1)
    
    fit_line_y = rate_eq(fit_line_x, k=rate, a=asymptote)
    make_error_lines = True
    try:
        fit_line_y_minus_error = rate_eq(fit_line_x, k=rate - confint, a=asymptote)
        fit_line_y_plus_error = rate_eq(fit_line_x, k=rate + confint, a=asymptote)
    except:
        make_error_lines = False
    
    main_line_symbol = 'k-'
    data_points_symbol = 'ro'
    plt.plot(fit_line_x, fit_line_y, main_line_symbol)
    
    if make_error_lines and not np.isnan(fit_line_y_plus_error).any() and not np.isnan(fit_line_y_minus_error).any():
        fit_line_y_plus_error[fit_line_y_plus_error > MAXIMUM_GRAPH_RATE_ERROR] = MAXIMUM_GRAPH_RATE_ERROR
        fit_line_y_minus_error[fit_line_y_minus_error < MINIMUM_GRAPH_RATE_ERROR] = MINIMUM_GRAPH_RATE_ERROR
    
    plt.plot(x, y, data_points_symbol)
    
    if make_error_lines:
        plt.fill_between(fit_line_x, fit_line_y_minus_error, fit_line_y_plus_error, color='black', alpha=.15)
    
    plt.tight_layout()
    
    plt.savefig("D:\\DR Testing\\graph.png", format="png")


if __name__ == "__main__":
    main()
    