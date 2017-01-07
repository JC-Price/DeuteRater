"""
Copyright (c) 2016 Bradley Naylor, Michael Porter, J.C. Price, and Brigham Young University

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

import numpy as np
import pandas as pd
import os
def normalize(data):
    norm_data = []
    total = sum(data)
    for x in data:
        norm_data.append(x/total)
    return norm_data


#does a check and removes outliers from neutromer spacing or combined
def remove_outlier(turnovers, cutoff, filters):
    med = np.median(turnovers)
    diff = []
    for point in  turnovers:
        diff.append(abs(point-med))
    mad = np.median(diff)
    zScore = []
    for item in diff:
        z = .6745 *item/mad
        zScore.append(z)
    true_turnovers =[]
    for pos in range(len(zScore)):
        if zScore[pos] > filters["ZCUTOFF"]:
            turnovers[pos] = "point {0} is an outlier and was removed".format(turnovers[pos])
        else:
            true_turnovers.append(turnovers[pos])
    if len(true_turnovers) > 1:
        return np.average(true_turnovers), np.std(true_turnovers, ddof =1), np.median(true_turnovers), turnovers 
    else:
        return np.average(turnovers), np.std(turnovers, ddof =1), np.median(turnovers), turnovers 


def correct_header(header, settings, filters):
    header.append('Experimental_Amount_of_Label (eaol) (0 = {0})'.format(filters["Zero Labeling Check"]))
    size = int(filters["Peaks included if over mass cutoff"])
    if settings["Use Abundance"]: 
        for i in range(size): header.append("M{0} difference from no labeling theory".format(i))
        for i in range(size): header.append("M{0} fraction new".format(i))
        header.extend(['average abundance fraction new','abundance std_dev', 'M0 fraction new for calculation'])
    if settings["Use neutromer spacing"]:
        for i in range(1, size): header.append("M{0}-M0 difference from no labeling theory".format(i))
        for i in range(1, size): header.append("M{0}-M0 fraction new".format(i))
        header.extend(['average neutromer spacing fraction new','neutromer spacing std_dev'])
        header.append("Points after outlier exclusion".format(i))
        header.append("median neutromer spacing fraction new")
    if settings["Use Abundance"] and settings["Use neutromer spacing"]:header.extend(['average combined fraction new','combined std_dev', 'median combined fraction new'])
    return header
    
    
#does the actual calculation of change in the 
def mida_D2O_calc(data, water, settings, filters):
    #set sizes
    #pandas series hold numpy floats, useful and flexible, but need to be ints here
    size = int(filters["Peaks included if over mass cutoff"])
    small_size = int(filters["Peaks included if under mass cutoff"])
    if data["theoretical mass"] > filters["mass cutoff"]:
        use_size = size
    else: use_size = small_size
    #put the amount of water into the main dataframe
    data['Experimental_Amount_of_Label (eaol) (0 = {0})'.format(filters["Zero Labeling Check"])]  = water['labeling'][data['file']]
    total_points = [] #if both abundance and neutromer spacing are used we need to track all points.
    if settings["Use Abundance"]:
        #makes the mida object
        row_mida = data['abundances'].split()
        row_mida = [float(s) for s in row_mida[:use_size]]
        normalized = normalize(row_mida)
        fraction_new = []
        for i in range(len(normalized)):
            max_delta = (data["M{0} maximum theory".format(i)]-data["M{0} no labeling theory".format(i)])
            diff =normalized[i] - data["M{0} no labeling theory".format(i)]
            data['M{0} difference from no labeling theory'.format(i)] = diff
            change = diff/max_delta
            if abs(max_delta) <= filters["Minimum Abund Change"]:
                data["M{0} fraction new".format(i)] = "Value {0} was not used due to low maximum possible change".format(change)
            else:
                data["M{0} fraction new".format(i)] = change
                fraction_new.append(change)
        total_points.extend(fraction_new)
        if use_size == small_size:
            for r in range(small_size, size):
                data['M{0} fraction new'.format(r)] = ''
        #there is a cutoff in the calculation.  the program will throw a fit if statistics are tried on empty lists so this controls that
        if fraction_new == []:
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] = 'insufficient labeling', 'insufficient labeling', 'insufficient labeling'
        elif type(data["M0 fraction new"]) == str:
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] ="M0 could not be used", "M0 could not be used", "M0 could not be used"
        elif len(fraction_new) == 1:  
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] =fraction_new[0], 'insufficient labeling', fraction_new[0]
        else:
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] =np.average(fraction_new), np.std(fraction_new,ddof =1), fraction_new[0]
    if settings["Use neutromer spacing"]:
        #similar to above analysis
        row_mida =data['observed m/z'].split()
        row_mida = [float(x) for x in row_mida[:use_size]]
        spacing_points = []
        for i in range(1,len(row_mida)):
            exp_diff = row_mida[i]-row_mida[0]
            max_delta = (data["M{0}-M0 maximum theory".format(i)] - data["M{0}-M0 no labeling theory".format(i)])
            delta = exp_diff - data["M{0}-M0 no labeling theory".format(i)]
            data["M{0}-M0 difference from no labeling theory".format(i)] = delta
            fraction_new = delta/max_delta
            data["M{0}-M0 fraction new".format(i)] = fraction_new
            spacing_points.append(fraction_new)
        total_points.extend(spacing_points)
        if use_size == small_size:
            for r in range(small_size,size):
                data["M{0}-M0 fraction new".format(r)] = ''
        if spacing_points == []:
            data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] = 'insufficient points', 'insufficient points', 'insufficient points'
            data['Points after outlier exclusion'] = "insufficient points"
        elif len(spacing_points) ==1 or len(spacing_points) ==2:
            if len(spacing_points) == 1:
                data['neutromer spacing std_dev'] = 'insufficient points'
                data['average neutromer spacing fraction new'], data['median neutromer spacing fraction new'] = spacing_points[0], spacing_points[0]
            else:
                data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] =np.average(spacing_points), np.std(spacing_points,ddof =1), np.median(spacing_points)
            data["Points after outlier exclusion"] = "Cannot perform outlier check on less than 3 points"
        else:
            #need to remove outliers (we trust the abundance is good above a certain cutoff
            ave, std, med, points = remove_outlier(spacing_points, filters["neutromer spacing Agreement Filter"], filters)
            data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] = ave, std, med
            data["Points after outlier exclusion"] = points
    if settings["Use Abundance"] and settings["Use neutromer spacing"]:
        #remove all outliers for this case(though mainly neutromer spacing)
        a,b,c, d= remove_outlier(total_points, filters["Combined Agreement Filter"], filters)
        if len(d) > 0:
            data["average combined fraction new"], data["combined std_dev"], data["median combined fraction new"], data["combined points"] = a, b, c, d
        else:
            data["average combined fractin new"], data["combined std_dev"], data["median combined fraction new"], data["combined points"] = np.nan, np.nan, np.nan, np.nan    
    return data



def theory_prep(theoretical, settings, size):
    #generate theoretical fitting output
    theory_header = ["sequence", "Elemental Composition", "theoretical mass", 'number of possible labeling sites', 'labeling']
    if settings["Use Abundance"]:
        for r in range(size):
            theory_header.append("M{0} no labeling theory".format(r))
        for r in range(size):
            theory_header.append("M{0} maximum theory".format(r))
    if settings["Use neutromer spacing"]:
        for r in range(1, size):
            theory_header.append("M{0}-M0 no labeling theory".format(r))
        for r in range(1, size):
            theory_header.append("M{0}-M0 maximum theory".format(r))
    Theory =  pd.DataFrame(list(theoretical), columns = theory_header)    
    Theory.to_csv(os.path.join(settings["Output Folder"],'Theoretical_Curves_for_Peptides.csv'), index = False)
    return Theory 
            
def Calculate(theoretical, experimental, water, settings, filters):
    #theory needs a header, transformation into a DataFrame, to be put in a file, and trimmed.  this accomplishes that.
    Theory = theory_prep(theoretical, settings, int(filters["Peaks included if over mass cutoff"]))
    #combine theory and experimental data together
    Complete = experimental.merge(Theory, on = ['sequence', 'labeling'])
    #filters out the peaks with insufficient labeling
    Complete = Complete[Complete['number of possible labeling sites'] >= filters["Required Number of Labels"]]
    #charge is not considered in the calculation of the theoretical curves (for abundance it is irrelevant and for neutromer spacing there may be multiple curves)
    #the correction is below: divide the initial values and all coefficients by the charge.  This has not shown issues.
    if settings["Use neutromer spacing"]:
        Complete.ix[:,'M1-M0 no labeling theory':] = Complete.ix[:,'M1-M0 no labeling theory':].astype(float).div(Complete['charge'].astype(float), axis ='index')
    Complete.to_csv(os.path.join(settings["Output Folder"],'Combination_of_Experimental_Data_and_Theory.csv'), index = False)
    header = list(Complete.columns.values) #for some reason the apply alphabetizes the columns, which is unhelpful.  we will save the header, update it and force the columns into line
    header = correct_header(header, settings, filters)
    #run the calculation of % labeled on the previous data and print it out
    Complete= Complete.reindex(columns = header)
    Complete = Complete.apply(mida_D2O_calc, axis= 1, args = (water, settings, filters))
    Complete.to_csv(os.path.join(settings["Output Folder"],'Calculation_of_Fraction_New_Protein.csv'), columns = header, index = False)