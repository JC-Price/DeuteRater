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

translated for theoretical kinetic calculation by Bradley Naylor with the help
of Spencer Lofthouse as of 11/12/14

updated to more current software development practices by Kyle Cutler and
Russell Denton as of 15 March 2019
"""

# NOTE: lines 200 and 201 need to be modified for each molecule..
# NOTE: could be sped up with c/c++ bindings to the original

from collections import namedtuple
from copy import deepcopy
import re


# name comes from Super Atom Data in the original emass program,
# this is the master list that holds all data throughout the calculations
sad = []

# electrons only matter for charge.
# Keep in mind this is for normal chemistry(charge means electrons move)
# ESI mass spectrometry may need changes
ELECTRON_MASS = .00054858
# obvious error mass to indicate issues
DUMMY_MASS = -10000000


# takes a sequence string(C2H5 for example) and turns it into a dictionary
# ( {'C':2, 'H':5} )
#does allow multiple letter symbols like Cu
#does NOT check the number of letters or if the letters are a valid chemical
#$does NOT allow leaving out 1 so if there is one sulfur 
#$that must be S1 not S
def new_parser(input_string):
    formmap = {}
    elements_separated = re.findall('[A-Z][^A-Z]*',input_string) 
    for e in elements_separated:
        element = e.rstrip('0123456789')
        number = e[len(element):]
        if number == "": number =1 #$ need the one to be there if element is alone
        formmap[element] = int(number)
    return formmap

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


# f is a pattern (list of peaks) limit is what a relative value must be under
# to be kept in.
#
# default is zero if you want it higher, change limit(near the emass function)
def prune(f, limit):
    prune = []
    counter = 0
    for i in f:
        if i.abundance > limit:
            break
        prune.append(counter)
        counter += 1
    counter = len(f) - 1
    for i in reversed(f):
        if i.abundance > limit:
            break
        prune.append(counter)
        counter -= 1
    for i in reversed(sorted(list(set(prune)))):
        del f[i]


# print_pattern is the output pattern.  name is a holdover from original emass,
# it does not print anymore
#
# the function normalizes the results so the highest is 100 (change if
# possible, extra function is inefficient)
#
# it will only add a value if the normalized value is higher than print_limit
# which is 5 *10^-7 by default
def print_pattern(result, digits):
    max_area = 0
    sum_area = 0
    for i in result:
        if max_area < i.abundance:
            max_area = i.abundance
        sum_area += i.abundance
    if max_area == 0:
        return
    print_limit = pow(10.0, -digits) / 2
    mass_list = []
    amount_list = []
    for i in result:
        val_perc = i.abundance / max_area * 100
        if i.mass != DUMMY_MASS and val_perc >= print_limit:
            mass_list.append(i.mass)
            amount_list.append(val_perc)
    return mass_list, amount_list


# this function calculates the isotope pattern by using convolute basic a
# bunch of times.
# tmp is an empty list for calculations(swapped with result when added to).fm
def calculate(result, fm, limit, charge):
    for i in fm:
        sal = [deepcopy(master_isotope[i])]
        n = int(fm[i])
        j = 0

        # this while loop is run on the number of atoms of a certain type.
        # reduced by a rightward bit shift.
        while n > 0:
            sz = len(sal)
            if j == sz:
                sal.append([])
                sal[j] = convolute_basic(sal[j - 1], sal[j - 1])
                prune(sal[j], limit)
            if n & 1:
                tmp = convolute_basic(result, sal[j])
                prune(tmp, limit)
                # don't want to copy all elements, requires fixing
                tmp, result = result, tmp
            n = n >> 1
            j += 1
    # this code must be uncommented to use the charge state.
    # adjustments are made
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

# X starts at normal isotope frequencies(in this case H.
# X is positions that can be deuterated, not will be TODO:??
#$ D is a label applied artificially (as in a standard) and does not change
master_isotope = {
    'X': [isotope(mass=1.0078246, abundance=0.999844),
          isotope(mass=2.0141021, abundance=0.000156)],
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
    # CQ Added to account for Flourine occurring in Lipids
    'F': [isotope(mass=18.99840322, abundance=1.0)],
    'D': [isotope(mass = 2.0141021, abundance=0.000156)],
    'Cl':[isotope(mass= 34.9688527, abundance=0.7576),
          isotope(mass=-1000000, abundance=0),
          isotope(mass= 36.9659026, abundance=0.2424)],
    'Br':[isotope(mass= 78.9183376, abundance=0.5069),
          isotope(mass=-1000000, abundance=0),
          isotope(mass= 80.9162897, abundance=0.4931)],
    'I': [isotope(mass=126.9044719, abundance=1.0)],
    'Si':[isotope(mass=27.9769265, abundance=0.92223),
          isotope(mass=28.9764946, abundance=0.04685),
          isotope(mass=29.9737701, abundance=0.03092)]
}

# TODO: What are these values and why are they here?
limit = 0
digits = 6
charge = 0

# will run emass return a list of relative abundances (4 or 5 things long) and
# a list of m/z(same peaks)
#
# note that these m/z values are essentially 'neutral' masses having their
# endogenous charge.
#
# the protons need to be added and the mass divided by charge to get m/z
# Phosphotidyl Serine elemental composition C46H77N1O10P1, mw = 834.stuff


def emass(cf, num_peaks,outfile=None):
    """
    Parameters
    ----------
    cf : str
        The chemical formula of the analyte
    num_peaks : int
        The number of expected isotope peaks
    Returns
    -------
    intensity list for the for num_peaks peaks with no enrichment
    """
    trunc_len = num_peaks  # This variable is for truncating lists.
    
    formatted_cf = new_parser(cf)
    result = calculate([isotope(0, 1)], formatted_cf, limit, charge)
    mz_list, intensity_list = print_pattern(result, digits)
    intensity_list = normalize(intensity_list[:trunc_len])
        
    return intensity_list




def main():

    intensity_values, true_m0 = emass(
        "C46H79N1O10P1", 5)
    print (intensity_values)
    print (true_m0)
    print("This is an auxiliary file to the n-value calculator, and should be treated as a library module")  # noqa


    return


if __name__ == "__main__":
    main()
