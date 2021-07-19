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


def calc_cf(aa_counts, aa_composition_df):
    elem_dict = {elem: 0 for elem in aa_composition_df.loc[:, 'C':].columns}
    for aa, count in aa_counts.items():
        for elem in elem_dict.keys():
            elem_dict[elem] += count * aa_composition_df.at[aa, elem]
    elem_dict = positive_cf_change(elem_dict)
    return elem_dict

#the calc_cf function calculates elemental composition well
#however it assumes all AAs have pepitde bonds on both ends  
#(which is not true of terminal AAs) and does not include the charge values
#so h2O and one H per charge must be added to get the most accurate cf
#if we go to negative mode or use different adducts, adjust here

#adds 1 O and z+2 H values to the elemental dict.
#takes elemental_dict (dictionary) and z (an int)
#returns the updated elemental dictionary
def positive_cf_change(elemental_dict):
    # since we are assuming positve and +H adduct we can just be direct
    elemental_dict["O"] +=1
    elemental_dict["H"] += 2# + z
    return elemental_dict


def calc_theory_mass(cf, elements_df):
    mass = 0
    for elem, count in cf.items():
        mass += count * elements_df.at[elem, 'relative_atomic_mass']
    return mass


def calc_add_n(aa_counts, aa_labeling_dict):
    add_n = 0
    for aa, count in aa_counts.items():
        add_n += count * aa_labeling_dict[aa]
    return add_n
