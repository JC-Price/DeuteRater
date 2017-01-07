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
import pandas as pd
import Isotopomer_Class_vrs2 as IC
import operator
import csv
import os

#get modifications
master_mod = {} #dictionary of the modifications
with open (os.path.join(os.path.dirname(__file__), "Modifications.csv"), "rt") as infile:
    reader = csv.reader(infile)
    next(reader, None)
    for row in reader:
        master_mod[row[0]] = row[1]
#opens a csv and pulls only the columns that are requested.  returns the header and all the data
def selective_readfile(location, columns):
    data = []
    with open(location, 'rt') as infile:
        reader = csv.reader(infile)
        for row in reader:
            selected_data = []
            try:
                for pos in columns:
                    selected_data.append(row[pos])
                data.append(selected_data)
            #ignores the row if the data is not present.
            except IndexError:
                y = 'insufficient data'
        #header and then data
        return data[0], data[1:]

#analyzes the modifications
def analyze_mods(mod, seq):
    #split the modification cell into its individual modifications
    mods = mod.split(';')
    if len(mods[0]) == 0:
        return seq
    else:
        temp_dict = {}
        for r in mods:
            #everything after the @ is the computers index
            loc = r.index('@')
            try:
                try: 
                    AA_key = int(r[loc+1:])-1
                except ValueError:#string for location check
                    if r[loc+1:] == "N-term" :
                        AA_key = 0
                    else: return ''
                temp_dict[AA_key] = master_mod[r[:loc]]#subtract 1 to get to the computer indexing
            except KeyError: # catch
                return '' # if the term is unknonw inerpretation is not possible so return the empty string to say dump the line
        try:
            newSeq = ''
            for i in range(len(seq)):
                if i in temp_dict.keys():
                    newSeq = newSeq + temp_dict[i]
                else:
                    newSeq = newSeq + seq[i]
        #bad location in sequence, dump the sequence
        except IndexError: 
            newSeq = ''
        return newSeq
        
#this function pulls the Protein ID out and reorders the list  data into the form that will be used for the rest of the program
def clean_up_data(data, f,AA_DICT, settings, filters):
    for r in reversed(range(len(data))):
        data[r][4] = analyze_mods(data[r][7],data[r][4])#must do first or the weight agreement will kill the peptide
        #this checks that there are enough peaks and that the thoeretical and measured m/z match. if not something is wrong and the row is deleted
        try:
            theory_wt = IC.peptide(data[r][4], settings["Heavy Element"], AA_DICT).MW()
            if theory_wt > filters["mass cutoff"]: target = filters["Peaks included if over mass cutoff"]
            else:target = filters["Peaks included if under mass cutoff"]
            x = [float(y) for y in data[r][-2].split()]
            x = [a for a in x if a != 0] # returns values as zero if not present.  remove issues.
            if len(x) < target or abs(theory_wt-float(data[r][6])) > filters["Weight Agreement"]or data[r][4] == '':#deals with insufficient peaks, bad weight agreement, and bad mods
                del data[r]
            else:
                #Protein ID, file, charge, RT, Peptide, protein name, mass, mod, isotopes, abundances before this point
                data[r] = [f, data[r][0], data[r][5], data[r][4], data[r][6], data[r][2], data[r][3], data[r][-3], data[r][-2], data[r][-1]] # if all is good set the data as good and proceed
        #bad sequences are removed
        except KeyError: 
            del data[r]


#data is the list, locations is a list of positions in the list, the filters should be floats, rt_ave is set for an agilent q-tof.  if you have a different machine please adjust
def initial_filter(data, locations, mz_filter, rt_filter, rt_ave = .5): #mz_filter = 15 ppm should replace this at some point
    #get the specific locations out so use more clearly
    pos_mz, pos_seq, pos_rt = tuple(locations)
    data_size = len(data)
    #need to be floats for the math to work
    for r in range(data_size):
        data[r][pos_mz] = float(data[r][pos_mz])
        data[r][pos_rt] = float(data[r][pos_rt])
    data = sorted(data, key=operator.itemgetter(pos_mz)) # this step allows us to filter more easily.  only have to check one way and until we are far enough away
    del_later = []

    for r in range(data_size):
        #start at the current positon and potentially go to the end of the list
        for i in range(r+1, data_size):
            #check if we are within m/z to start with to limit the number of things the computer must do
            if data[i][pos_mz] <= data[r][pos_mz] + mz_filter:
                if data[i][pos_seq] != data[r][pos_seq]:
                    if data[r][pos_rt] - rt_filter<= data[i][pos_rt] <= data[r][pos_rt] + rt_filter:
                       #this is needed because we don't know how many matches there are, so we just save every position(deleting in place is trickier)
                       del_later.append(r)
                       del_later.append(i)
            else:break
    #no replicates and sort it. then delete positions.
    del_later = sorted(list(set(del_later)))
    for position in reversed(del_later):
        del(data[position])
    return data

#this function filters out a particular issue.  it is possible that extraction may leave out peaks. 
def check_mz(data, mz_pos, charge_pos, abund_pos, mass_pos, filters):
    filtered_data = []
    for row in data:
        charge = float(row[charge_pos])
        mz = row[mz_pos].split(' ')
        abund = row[abund_pos].split(' ')
        for x in range(len(mz)):
            mz[x] = float(mz[x])
            abund[x] =float(abund[x])
        r =1 
        good = True
        #delete bad peaks 
        while r < len(mz):
            #each peak should be within 1 m/z unit (+ or - .1 m/z unit) to be real.
            #yes this is a wide and arbitrary window.  This is not designed to filter out bad matches, just get rid of data where the pattern is incomplete
            #full filtering should be done in the extractor 
            if mz[r] > mz[r-1] + 1/charge +.1 or mz[r] < mz[r-1] + 1/charge - .1:
                del mz[r]
                del abund[r]
                good = False
            else:
                r += 1
        if good:
            filtered_data.append(row)
        # this deals with the case of extra peaks added inbetween good peaks (will probably never happen)
        elif r >= filters["Peaks included if over mass cutoff"] or (float(row[mass_pos]) < filters["mass cutoff"] and r ==filters["Peaks included if under mass cutoff"]):
            mz_string = ''
            abund_string  = ''
            for r in range(len(mz)):
                mz_string += '{0} '.format(mz[r])
                abund_string += '{0} '.format(abund[r])
                row[mz_pos] = mz_string[:-1]
                row[abund_pos] = abund_string[:-1]     
            filtered_data.append(row)   
    return filtered_data     

def begin_preparation(settings, filters, allowed_files, AA_DICT):
    #make needed data structures
    experimental = {}
    panda = pd.DataFrame()
    #make compressed files directory if necessary
    if not os.path.isdir(os.path.join(settings["Output Folder"], "Compressed_Datafiles")):
        os.makedirs(os.path.join(settings["Output Folder"], "Compressed_Datafiles"))
    #cycles through all files in input directory and extracts data
    bad_files = ""
    for f in allowed_files:
        try:
            #get data
            header, data = selective_readfile(os.path.join(settings["Input Folder"], f), [0,1, 2,3,4,5,6, 8,12, 13,16])
            #need to reverse list to delete rows if need be
            #this particular check is to ensure that the m/z agrees with what the elemental composition says and that there are sufficient peaks for later analysis.
            clean_up_data(data,f, AA_DICT, settings, filters) #puts data into proper format for the rest of the program
            data = initial_filter(data, [4,3,6], filters["MZ Filter"], filters["RT Filter"])# filters sequences that are too close in both rt and m/z (too likely that a peak will be incorrectly assigned in calculations)
            data = check_mz(data, 7, 5, 8, 4, filters)#filters sequences where some peaks are missed.
            header = ['file', 'Protein ID', "Protein name", 'sequence', 'experimental mass', 'charge', 'observed retention time', 'observed m/z', 'abundances', "homologous proteins"]
            data = pd.DataFrame(data, columns = header)
            experimental[f] = data
            #puts the compressed and filtered data into a new file for the user to look at.
            data.to_csv(os.path.join(settings["Output Folder"], 'Compressed_Datafiles', '{}_compressed.csv'.format(f)), index = False)
            panda = panda.append(data)
        except:
            bad_files += f
            bad_files += ", "
    if len(bad_files) > 0:
        bad_files = bad_files[:-2]
    return panda, bad_files