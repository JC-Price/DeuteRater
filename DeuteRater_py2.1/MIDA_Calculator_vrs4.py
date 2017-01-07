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
import Mida_Calculation_Class_vrs2 as mcc
import os
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
        for i in range(size): header.append("M{0} delta maximum from eaol".format(i))
        for i in range(size): header.append("M{0} fraction new".format(i))
        header.extend(['average abundance fraction new','abundance std_dev', 'M0 fraction new for calculation'])
    if settings["Use neutromer spacing"]:
        for i in range(1, size): header.append("M{0}-M0 difference from no labeling theory".format(i))
        for i in range(1, size): header.append("M{0}-M0 delta maximum from eaol".format(i))
        for i in range(1, size): header.append("M{0}-M0 fraction new".format(i))
        header.extend(['average neutromer spacing fraction new','neutromer spacing std_dev'])
        header.append("Points after outlier exclusion")
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
        row_mida = mcc.mida(data['theoretical mass'], data['abundances'].split(), data['M0 coefficient a':'delta M{0} cubic coeff 1'.format(use_size-1)], data['Experimental_Amount_of_Label (eaol) (0 = {0})'.format(filters["Zero Labeling Check"])], data['M0 no labeling theory':'M{0} no labeling theory'.format(use_size-1)], True, filters)
        row_turnover, row_reporter = row_mida.change() # does the calculations
        total_points.extend(row_turnover)
        #makes a pandas series for easy addition into the main dataframe
        for r in range(len(row_mida.delta)):
            data['M{0} difference from no labeling theory'.format(r)] = row_mida.delta[r]
        if use_size == small_size:
            for r in range(small_size, size):
                data['M{0} difference from no labeling theory'.format(r)] = ''
        for r in range(len(row_mida.max_change)):
            data['M{0} delta maximum from eaol'.format(r)] = row_mida.max_change[r]
        #make blanks for formatting 
        if use_size == small_size:
            for r in range(small_size, size):
                data['M{0} delta maximum from eaol'.format(r)] = ''
        
        for r in range(len(row_reporter)):
            data['M{0} fraction new'.format(r)] = row_reporter[r]
        #make blanks for formatting 
        if use_size == small_size:
            for r in range(small_size, size):
                data['M{0} fraction new'.format(r)] = ''
        #there is a cutoff in the calculation.  the program will throw a fit if statistics are tried on empty lists so this controls that
        if row_turnover == []:
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] = 'insufficient labeling', 'insufficient labeling', 'insufficient labeling'
        elif len(row_turnover) == 1:  
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] =np.average(row_turnover[0]), 'insufficient labeling', row_turnover[0]
        elif type(row_reporter[0]) != str:
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] =np.average(row_turnover), np.std(row_turnover,ddof =1), row_turnover[0]
        else:
            value_loc = max(row_mida.max_change)
            loc = row_mida.max_change.index(value_loc)
            data['average abundance fraction new'], data['abundance std_dev'], data['M0 fraction new for calculation'] =np.average(row_turnover), np.std(row_turnover,ddof =1), row_reporter[loc] # in this case M0 is not usable, so row_turnover is small.  Need an index with everything (we have one or more things that are valid or we would be caught in the first if)
    if settings["Use neutromer spacing"]:
        #similar to above analysis
        row_mida = mcc.mida(data['theoretical mass'], data['observed m/z'].split(), data['M1-M0 4th order coeff 4': 'M{0}-M0 4th order coeff 1'.format(size-1)], data["Experimental_Amount_of_Label (eaol) (0 = {0})".format(filters["Zero Labeling Check"])], data['M1-M0 no labeling theory':'M{0}-M0 no labeling theory'.format(size-1)],False,filters)
        row_turnover, row_reporter = row_mida.change()
        total_points.extend(row_turnover)
        for r in range(len(row_mida.delta)):
            data['M{0}-M0 difference from no labeling theory'.format(r+1)] = row_mida.delta[r]   
        if use_size == small_size-1: 
            for r in range(small_size, size):
                data['M{0}-M0 difference from no labeling theory'.format(r)] = ''
        for r in range(len(row_mida.max_change)):
            data['M{0}-M0 delta maximum from eaol'.format(r+1)] = row_mida.max_change[r]
        if use_size == small_size-1: 
            for r in range(small_size,size):
                data['M{0}-M0 delta maximum from eaol'.format(r)]= ''
        for r in range(len(row_reporter)):
            data['M{0}-M0 fraction new'.format(r+1)] = row_reporter[r]
        if use_size == small_size-1: 
            for r in range(small_size,size):
                data['M{0}-M0 fraction new'.format(r)]= ''
        if row_turnover == []:
            data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] = 'insufficient points', 'insufficient points', 'insufficient points'
            data['Points after outlier exclusion'] = "insufficient points"
        elif len(row_turnover) == 1 or len(row_turnover) ==2:
            if len(row_turnover) == 1:
                data['neutromer spacing std_dev'] = 'insufficient points'
                data['average neutromer spacing fraction new'], data['median neutromer spacing fraction new'] = row_turnover[0], row_turnover[0]
            else:
                data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] =np.average(row_turnover), np.std(row_turnover,ddof =1), np.median(row_turnover)
            data["Points after outlier exclusion"] = "Cannot perform outlier check on less than 3 points"
        else:        
            #need to remove outliers (we trust the abundance is good above a certain cutoff
            ave, std, med, points = remove_outlier(row_turnover, filters["neutromer spacing Agreement Filter"], filters)
            data['average neutromer spacing fraction new'], data['neutromer spacing std_dev'], data['median neutromer spacing fraction new'] = ave, std, med
            data['Points after outlier exclusion'] = points
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
    theory_header = ["sequence", "Elemental Composition", "theoretical mass", 'number of possible labeling sites']
    if settings["Use Abundance"]:
        for r in range(size):
            theory_header.append("M{0} no labeling theory".format(r))
        theory_header.extend(["M0 coefficient a", "M0 coefficient k", "delta M0 pearson r"])
        for r in range(1, size):
            theory_header.extend(['delta M{0} cubic coeff 3'.format(r), 'delta M{0} cubic coeff 2'.format(r), 'delta M{0} cubic coeff 1'.format(r), 'delta M{0} pearson r'.format(r)])
    if settings["Use neutromer spacing"]:
        for r in range(1, size):
            theory_header.append("M{0}-M0 no labeling theory".format(r))
        for r in range(1, size):
            theory_header.extend(['M{0}-M0 4th order coeff 4'.format(r), 'M{0}-M0 4th order coeff 3'.format(r), 'M{0}-M0 4th order coeff 2'.format(r), 'M{0}-M0 4th order coeff 1'.format(r), 'M{0}-M0 pearson r'.format(r)])
    Theory =  pd.DataFrame(list(theoretical), columns = theory_header)    
    Theory.to_csv(os.path.join(settings["Output Folder"],'Theoretical_Curves_for_Peptides.csv'), index = False)
    #r values are not necessary for the rest of the calculations, and so they are removed after they are put in a csv
    drop_it =[]
    for word in Theory.columns:
        if word[-1] == 'r': drop_it.append(word)
    Theory = Theory.drop(drop_it, 1)
    return Theory 
            
def Calculate(theoretical, experimental, water, settings, filters):
    #theory needs a header, transformation into a DataFrame, to be put in a file, and trimmed.  this accomplishes that.
    Theory = theory_prep(theoretical, settings, int(filters["Peaks included if over mass cutoff"]))
    #combine theory and experimental data together
    Complete = experimental.merge(Theory, on = 'sequence')
    #filters out the peaks with insufficient labeling
    Complete = Complete[Complete['number of possible labeling sites'] >= filters["Required Number of Labels"]]
    #charge is not considered in the calculation of the theoretical curves (for abundance it is irrelevant and for neutromer spacing there may be multiple curves)
    #the correction is below: divide the initial values and all coefficients by the charge.  This has not shown issues.
    if settings["Use neutromer spacing"]:
        Complete.ix[:,'M1-M0 no labeling theory':] = Complete.ix[:,'M1-M0 no labeling theory':].astype(float).div(Complete['charge'].astype(float), axis ='index')
    Complete.to_csv(os.path.join(settings["Output Folder"],'Combination_of_Experimental_Data_and_Theoretical_Curves.csv'), index = False)
    header = list(Complete.columns.values) #for some reason the apply alphabetizes the columns, which is unhelpful.  we will save the header, update it and force the columns into line
    header = correct_header(header, settings, filters)
    #run the calculation of % labeled on the previous data and print it out
    Complete= Complete.reindex(columns = header)
    Complete = Complete.apply(mida_D2O_calc, axis= 1, args = (water, settings, filters))
    Complete.to_csv(os.path.join(settings["Output Folder"],'Calculation_of_Fraction_New_Protein.csv'), columns = header, index = False)
def Recalculate(settings, filters):
    try:
        water = pd.read_csv(os.path.join(settings["Output Folder"],"Time_and_Labeling_Data.csv"))
        try:
            water = water.set_index('file')
        except:
            return "Cannot find 'file' in header. Please type 'file' in the header of the first column and Try again"    
        Complete = pd.read_csv(os.path.join(settings["Output Folder"], 'Combination_of_Experimental_Data_and_Theoretical_Curves.csv'))
    except IOError:
        return "Recalculation could not be completed, correct files were not present.  Check that Time_and_Labeling_Data.csv and Combination_of_Experimental_Data_and_Theoretical_Curves.csv are present in your Analysis Output Folder."
    try:
        a = water['labeling']
    except:
        return "Cannot find 'labeling' in header. Please type 'labeling' in proper column and Try again'"
    if water["labeling"].dtype != np.float64 and water["labeling"].dtype != np.int64:
        return "Non-numberical value for labeling in Time_and_Labeling_Data.csv.  Correct and Try again"
    for value in water["labeling"]:
        if np.isnan(value):
            return "Non-existant number for labeling in Time_and_Labeling_Data.csv.  Correct and Try again"
        if value > 1:
            return "Labeling value above 1 in Time_and_Labeling_Data.csv.  Correct and Try again"
        if value < 0:
            return "Labeling value below 0 in Time_and_Labeling_Data.csv.  Correct and Try again"        
    try:
        header = list(Complete.columns.values) #for some reason the apply alphabetizes the columns, which is unhelpful.  we will save the header, update it and force the columns into line
        if "M0 coefficient a" in header:
            header = correct_header(header, settings, filters)
            #run the calculation of % labeled on the previous data and print it out
            Complete= Complete.reindex(columns = header)
            Complete = Complete.apply(mida_D2O_calc, axis= 1, args = (water, settings, filters))
            Complete.to_csv(os.path.join(settings["Output Folder"], 'Calculation_of_Fraction_New_Protein.csv'), cols = header, index = False)
            return "Recalculation completed successfully"
        else: 
            return "Recalculation could not be completed.  Fast calculation was used or header has been changed"
    except:
       return "Unknown Error occured.  Please check your settings and that all files are present in you Analysis Output Folder without alteration."