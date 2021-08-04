"""
translation of emass which was originally done in C++
which contained the following copyright and disclaimer:
Based on an algorithm developed by Alan L. Rockwood.
Published in
Rockwood, A.L. and Haimi, P.: "Efficient calculation of
Accurate Masses of Isotopic Peaks",
Journal of The American Society for Mass Spectrometry
JASMS 03-2263, in press
Copyright (c) 2005 Perttu Haimi and Alan L. Rockwood
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
NOTE: lines 200 and 201 need to be modified for each molecule..
translated for theoretical kinetic calculation by Bradley Naylor with the help of Spencer Lofthouse
as of 11/12/14
updated to more current software development practices by Kyle Cutler and Russell Denton as of 15 March 2019
"""

from collections import namedtuple
from copy import deepcopy
import pandas as pd
import numpy as np
import re

# name comes from Super Atom Data in the original emass program,
# this is the master list that holds all data throughout the calculations
sad = []

# electrons only matter for charge.
# Keep in mind this is for normal chemistry(charge means electrons move) ESI mass spectrometry may need changes
ELECTRON_MASS = .00054858
# obvious error mass to indicate issues
DUMMY_MASS = -10000000

#$allows use of multiple letters for elemental composition
def new_parser(input_string):
    formmap = {}
    elements_separated = re.findall('[A-Z][^A-Z]*',input_string) 
    for e in elements_separated:
        element = e.rstrip('0123456789')
        number = e[len(element):]
        if number == "": number =1 #$ need the one to be there if element is alone
        formmap[element] = int(number)
    return formmap
"""
# takes a sequence string(C2H5 for example) and turns it into a dictionary ( {'C':2, 'H':5} )
def parser(elem_comp):
    #NOTE: cannot handle multiple letter symbols like Al or Se.  Modify if such is needed
    # initialize values
    letter = ''
    num = -100
    numstr = ''
    formmap = {}  # formula map.  holdover from original program.
    for e in elem_comp:
        # if e is a letter we are done with numbers for the previous 3 letter
        if e.isalpha():
            # no point putting '':-100 into the dictionary
            if letter != '' and num >= 0:
                formmap[letter] = num
                num = -100
            letter = e
        elif str(e).isnumeric():  # str was unicode (python 3 changed unicode to str)
            if num < 0:
                numstr = e
            else:
                numstr += e
            num = float(numstr)
    # won't put last value in without this last command
    formmap[letter] = num
    return formmap
"""

# combines two patterns (lists of peaks).
# The idea is that small groups of atoms are combined into large atoms.
# This function does the combination
def convolute_basic(g, f):
    """NOTE: math poorly understood do not change"""
    h = []
    g_n = len(g)
    f_n = len(f)
    if g_n == 0 or f_n == 0:
        return
    for k in range(g_n + f_n - 1):
        sumweight = 0
        summass = 0
        start = 0 if k < f_n - 1 else k - f_n + 1
        end = k if k < g_n - 1 else g_n - 1
        for i in range(start, end + 1):
            weight = g[i].abundance * f[k - i].abundance
            mass = g[i].mass + f[k - i].mass
            sumweight += weight
            summass += weight * mass
        if sumweight == 0:
            p = isotope(DUMMY_MASS, sumweight)
        else:
            p = isotope(summass / sumweight, sumweight)
        h.append(p)
    return h


# f is a pattern (list of peaks) limit is what a relative value must be under to be kept in.
# default is zero if you want it higher, change limit(near the emass function)
def prune(f, limit):
    prune = []
    counter = 0
    for i in f:
        if i.abundance > limit: break
        prune.append(counter)
        counter += 1
    counter = len(f) - 1
    for i in reversed(f):
        if i.abundance > limit: break
        prune.append(counter)
        counter -= 1
    for i in reversed(sorted(list(set(prune)))):
        del f[i]


# print_pattern is the output pattern.  name is a holdover from original emass, it does not print anymore
# the function normalizes the results so the highest is 100 (change if possible, extra function is inefficient)
# it will only add a value if the normalized value is higher than print_limit which is 5 *10^-7 by default
def print_pattern(result, digits):
    max_area = 0
    sum_area = 0
    for i in result:
        if max_area < i.abundance:
            max_area = i.abundance
        sum_area += i.abundance
    if max_area == 0: return
    print_limit = pow(10.0, -digits) / 2
    mass_list = []
    amount_list = []
    for i in result:
        val_perc = i.abundance / max_area * 100
        if i.mass != DUMMY_MASS and val_perc >= print_limit:
            mass_list.append(i.mass)
            amount_list.append(val_perc)
    return mass_list, amount_list


# this function calculates the isotope pattern by using convolute basic a bunch of times.
# tmp is an empty list for calculations(swapped with result when added to).fm
def calculate(result, fm, limit, charge):
    for i in fm:
        sal = [deepcopy(master_isotope[i])]
        n = int(fm[i])
        j = 0
        # this while loop is run on the number of atoms of a certain type.  reduced by a rightward bit shift.
        while n > 0:
            sz = len(sal)
            if j == sz:
                sal.append([])
                sal[j] = convolute_basic(sal[j - 1], sal[j - 1])
                prune(sal[j], limit)
            if n & 1:
                tmp = convolute_basic(result, sal[j])
                prune(tmp, limit)
                tmp, result = result, tmp  # don't want to copy all elements, requires fixing
            n = n >> 1
            j += 1
    # this code must be uncommented to use the charge state adjustments are made
    # for i in result:
    #     if charge > 0:
    #         i.mass = i.mass/abs(charge) - ELECTRON_MASS
    #     elif charge < 0:
    #         i.mass = i.mass/abs(charge) + ELECTRON_MASS
    return result


# standard function to normalize the peaks by ratio
# (sum of values will be one, as opposed to normalizing to the maximum value).
def normalize(old_list):
    s = sum(old_list)
    return [float(val)/s for val in old_list]

def normalizeM0(old_list):
    M0 = old_list[0]
    return [float(val)/M0 for val in old_list]


# stores isotope data
isotope = namedtuple('isotope', 'mass abundance')

# X starts at normal isotope frequencies(in this case H. X is positions that can be deuterated, not will be TODO:??
#$ D is a label applied artificially (as in a standard) and does not change
master_isotope = {'X': [isotope(mass=1.0078246, abundance=0.999844),
                        isotope(mass=2.0141021, abundance=0.000156)],  # TODO: this is identical to hydrogen?
                  'H': [isotope(mass=1.0078246, abundance=0.999844),
                        isotope(mass=2.0141021, abundance=0.000156)],
                  'C': [isotope(mass=12.0000000, abundance=0.9891),
                        isotope(mass=13.0033554, abundance=0.0109)],
                  'N': [isotope(mass=14.0030732, abundance=0.99635),
                        isotope(mass=15.0001088, abundance=0.00365)],
                  'O': [isotope(mass=15.9949141, abundance=0.99759),
                        isotope(mass=16.9991322, abundance=0.00037),
                        isotope(mass=17.9991616, abundance=0.00204)],
                  'P': [isotope(mass=30.973762, abundance=1.0)],
                  'S': [isotope(mass=31.972070, abundance=0.9493),
                        isotope(mass=32.971456, abundance=0.0076),
                        isotope(mass=33.967866, abundance=0.0429),
                        isotope(mass=-1000000, abundance=0),
                        isotope(mass=35.967080, abundance=0.0002)],
                  'F': [isotope(mass=18.99840322, abundance=1.0)],
                  'D': [isotope(mass = 2.0141021, abundance=1.0)],
                  'Cl':[isotope(mass= 34.9688527, abundance=0.7576),
                        isotope(mass=-1000000, abundance=0),
                        isotope(mass= 36.9659026, abundance=0.2424)],
                  'Br':[isotope(mass= 78.9183376, abundance=0.5069),
                        isotope(mass=-1000000, abundance=0),
                        isotope(mass= 80.9162897, abundance=0.4931)],
                  'I': [isotope(mass=126.9044719, abundance=1.0)],
                  'Si':[isotope(mass=27.9769265, abundance=0.92223),
                        isotope(mass=28.9764946, abundance=0.04685),
                        isotope(mass=29.9737701, abundance=0.03092)],
                  'Na':[isotope(mass=22.98976928, abundance=1.0)],
                  
                  } # CQ Added to account for Flourine occurring in Lipids

# TODO: What are these values and why are they here?
limit = 0
digits = 6
charge = 0

# will run emass return a list of relative abundances (4 or 5 things long) and a list of m/z(same peaks)
# note that these m/z values are essentially 'neutral' masses having their endogenous charge.
# the protons need to be added and the mass divided by charge to get m/z
# Phosphotidyl Serine elemental composition C46H77N1O10P1, mw = 834.stuff


def emass(chemical_formula, n_begin, n_end, low_pct, high_pct, num_peaks, step=1, testing=False, outfile=None):
    """
    Parameters
    ----------
    chemical_formula : str
        The chemical formula of the analyte
    n_begin : int
        The lower bound of the possible number of deuterium incorporation sites
    n_end : int
        The upper bound of the possible number of deuterium incorporation sites
    low_pct : float
        The lower deuterium enrichment value from the LGR
    high_pct : float
        The higher deuterium enrichment value from the LGR
    num_peaks : int
        The number of expected isotope peaks
    step : int
        The step between n-values when iterating
    Returns
    -------
    (pd.Dataframe, pd.Dataframe), (pd.Dataframe, pd.Dataframe)
        These pandas dataframes contain data from the lower and higher water enrichment values, respectively
    """
    trunc_len = num_peaks  # This variable is for truncating lists. TODO: discuss what this means

    def populate_dataframes(pct, testing=False):
        mz_lists = []
        intensity_lists = []

        # Added for testing purposes ~ Chad Quilling
        M0_intensity_lists = []
        full_mz_lists = []
        full_intensity_lists = []
        full_M0_intensity_lists = []

        shrunk_mz_lists = []
        shrunk_intensity_lists = []

        master_isotope['X'] = [isotope(master_isotope['H'][0].mass, master_isotope['H'][0].abundance - pct),
                               isotope(master_isotope['H'][1].mass, master_isotope['H'][1].abundance + pct)]

        # Fix master_isotope['X'] so the abundances are never negative or over 1.
        if master_isotope['X'][0].abundance < 0.0:
            master_isotope['X'][0] = isotope(master_isotope['H'][0].mass, 0.0)
        if master_isotope['X'][1].abundance > 1.0:
            master_isotope['X'][0] = isotope(master_isotope['H'][1].mass, 1.0)

        shrunk_num = 0
        for n in range(n_begin, n_end, step):
            n_h = n_end
            chem_format = new_parser(chemical_formula.format(n_h - n, n))
            result = calculate([isotope(0, 1)], chem_format, limit, charge)
            mz_list, intensity_list = print_pattern(result, digits)
            # the lengths of these lists are 11+ before truncating
            mz_lists.append([n] + mz_list[:trunc_len])
            intensity_lists.append([n] + normalize(intensity_list[:trunc_len]))

            # Added for testing purposes ~ Chad Quilling
            # M0_intensity_lists.append([n] + normalizeM0(intensity_list[:trunc_len]))
            full_mz_lists.append([n] + mz_list[:])
            full_intensity_lists.append([n] + normalize(intensity_list[:]))

            if shrunk_num == 0:
                shrunk_num = len(mz_list)

            shrunk_mz_lists.append([n] + mz_list[:shrunk_num])
            shrunk_intensity_lists.append([n] + normalize(intensity_list[:shrunk_num]))

            # full_M0_intensity_lists.append([n] + normalizeM0(intensity_list[:]))

        if not testing:
            # Use normilzation by number of peaks in n_D = 0
            # return (pd.DataFrame(data=shrunk_mz_lists, columns=['n_D'] + ['mz' + str(i) for i in range(shrunk_num)]),
            #         #         # CHOOSE IF YOU WANT TO OUTPUT NORMALIZATION BY SUM OR M0:
            #         #
            #         #         # SUM
            #                 pd.DataFrame(data=shrunk_intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(shrunk_num)])

            return (pd.DataFrame(data=mz_lists, columns=['n_D'] + ['mz' + str(i) for i in range(trunc_len)]),
            # CHOOSE IF YOU WANT TO OUTPUT NORMALIZATION BY SUM OR M0:

            # SUM
                            pd.DataFrame(data=intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)]))

            # M0 Uncomment from above as well
            # pd.DataFrame(data=M0_intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)]))

        else:
            return (pd.DataFrame(data=mz_lists, columns=['n_D'] + ['mz' + str(i) for i in range(trunc_len)]),
                    pd.DataFrame(data=intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)]),
                    pd.DataFrame(data=M0_intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)]),
                    pd.DataFrame(data=full_mz_lists, columns=['n_D'] + ['mz' + str(i) for i in range(trunc_len)]),
                    pd.DataFrame(data=full_intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)]),
                    pd.DataFrame(data=full_M0_intensity_lists, columns=['n_D'] + ['I' + str(i) for i in range(trunc_len)])
                    )

    # print("cf == ", chemical_formula)  # logging purposes

    unlabeled_dfs = populate_dataframes(low_pct, testing)
    labeled_dfs = populate_dataframes(high_pct, testing)
    return unlabeled_dfs, labeled_dfs


def main():
    # uncomment lines below to step through the execution in a debugging tool and/or print example to console
    unlabeled_dfs, labeled_dfs = emass("C46H{}BrN1O10Cl20P1X{}", 5, 79, 0, .053, 4, step=1)
    (df_unlabeled_mzs, df_unlabeled_intensities) = unlabeled_dfs
    (df_labeled_mzs, df_labeled_intensities) = labeled_dfs
    print(df_unlabeled_intensities.to_string())
    print("This is an auxiliary file to the n-value calculator, and should be treated as a library module")

    # uncomment lines below to get outputs: Normalized by M0 and not truncated
    # unlabeled_dfs, labeled_dfs = emass("C46H{}N1O10P1X{}", 5, 79, 0, .053, 4, step=1, testing=True)
    # (df_unlabeled_mzs, df_unlabeled_intensities, df_M0_unlabeled_intensities, df_full_unlabeled_mzs, df_full_unlabeled_intensities, df_full_unlabeled_M0_intensities) = unlabeled_dfs
    # (df_labeled_mzs, df_labeled_intensities, df_M0_labeled_intensities, df_full_labeled_mzs, df_full_labeled_intensities, df_full_labeled_M0_intensities) = labeled_dfs
    # print(df_unlabeled_intensities.to_string())
    # print("This is an auxiliary file to the n-value calculator, and should be treated as a library module")

    return

if __name__ == "__main__":
    main()
