"""
Copyright (c) 2016 Michael Porter, Bradley Naylor, J.C. Price, and Brigham Young University

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

from pyteomics import mzml
from Identification_v3 import *
from Peak import *
from Isotope import *
import time
import bisect
import os
import numpy as np
import pandas as pd
import warnings


def extract(id_arg,spec_arg,out_arg,n_arg,lz_arg,hz_arg,tW,ppmW):
    warnings.filterwarnings("ignore")            
                 
    idfile = id_arg
    spectrumfile = spec_arg
    out = out_arg
    
    NEUTRON = 1.008664916
    PROTON = 1.007276467
    N_ISOS = n_arg # how many peaks to look for
    lowZ = lz_arg
    highZ = hz_arg
    timeWindow = tW
    ppmWindow = ppmW

    ##############################################
    # The following functions are used for       #
    # calculating the angle between two vectors  #
    ##############################################
    # Taken from http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python retrieved 30 September 2015
    def unit_vector(vector):
        return vector/np.linalg.norm(vector)

    def angle_between(v1,v2):
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        angle = np.arccos(np.dot(v1_u,v2_u))
        if np.isnan(angle):
            if(v1_u == v2_u).all():
                return 0.0
            else:
                return np.pi
        return angle

    #################################################
    # The following function is used for condensing #
    # the chromatographic series of isotopes into a #
    # single isotopic envelope.  This function is   #
    # open to modification and optimization.        #
    #################################################
    def condense(key,isoList):
        mz = []
        rt = []
        intensity = []
        isoList = sorted(isoList,key=lambda sort: sort.getRT())
        for thing in isoList:
            if thing.getType()==key:
                mz.append(thing.getMZ())
                rt.append(thing.getRT())
                intensity.append(thing.getI())
                length = len(mz)
        if length==0:
            length = 1
            intensity.append(0)

        #median abs deviation
        #http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
        mzMedian = np.median(mz)
        mzDiff = []
        for thing in mz:
            mzDiff.append(abs(thing-mzMedian))
        mad = np.median(mzDiff)
        zScore = []
        for thing in mzDiff:
            modZ = 0.6745*thing/mad
            zScore.append(modZ)
        newMZ = []
        for i in range(0,len(mz)):
            if zScore[i]<=2.0:
                newMZ.append(mz[i])
        return Peak(np.mean(newMZ),np.median(intensity),np.median(rt))

    startTime = time.time()
    idList = [] # for storing the peptide IDs

    # Time to read in the ID file
    # Eventually should switch over to Pandas dataframe so
    # the column order is not so stringent.
    outfile = open(out,'w') # where I'm writing the envelopes
    outfile.write('\"Title\",\"ID\",\"Charge\",\"Retention Time\",\"Peptide Sequence\",\"Protein\",\"Peptide mass\",\"Precursor mass\",\"Modifications\",\"MSMS Info\",\"Base Peak m/z\",\"Envelope Size\",\"Isotopic Envelope m/z values\",\"Isotopic Envelope abundance values\",\"Min RT\",\"Max RT\",\"homologous proteins\"\n')
    try:
        #Now I read the results file into a list
        print 'Reading in ID file...'
        if idfile[-4:] == ".tsv":
            df_psm_all = pd.read_csv(idfile,index_col=False,header=0,delimiter='\t')
        elif idfile[-4:] == ".csv":
            df_psm_all = pd.read_csv(idfile,index_col=False,header=0)
        df_psm = df_psm_all.drop_duplicates(subset=['Sequence','Variable Modification(s)','Peptide Theoretical Mass'])
        df_psm['Variable Modification(s)'] = df_psm['Variable Modification(s)'].astype(str)
        #use the index as the drop command earlier will have made a range or similar unworkable (index does not shift to fill in gaps from a deletion in pandas
        for i in df_psm.index: 
            newMods = ""
            #the str command turned all empty variable mods into nan, and the format needs more than that to function so this is reasonable
            if len(df_psm['Variable Modification(s)'][i])>3:
                types = df_psm['Variable Modification(s)'][i].split("), ")
                for thing in types:
                    name = thing.split(" ")[0]
                    pos = thing.split("(")[1]
                    temp = pos.split(")")[0]
                    positions = temp.split(", ")
                    for p in positions:
                        newMods = newMods+name+"@"+p+";"
                newMods = newMods[0:-1]
            title = str(df_psm['Protein ID'][i])
            homology = str(df_psm["Homologous Proteins"][i])
            ident = str(df_psm['Data File Name'][i])
            rt = str(df_psm['Precursor Retention Time (sec)'][i]/60.0)
            pep = str(df_psm['Sequence'][i])
            prot = str(df_psm['Protein Name'][i])
            mh = float(df_psm['Peptide Theoretical Mass'][i]+PROTON)
            tmpMass = float(df_psm['Precursor m/z'][i])
            tmpCharge = float(df_psm['Identification Charge'][i])
            pepmass = mh-PROTON
            premass = tmpMass*tmpCharge-tmpCharge*PROTON
            mods = newMods
            envelopesize = 0
            mz = []
            intensity = []
            for z in range(lowZ,highZ+1):
                basepeak = (premass+float(z)*PROTON)/float(z)
                mspre = basepeak
                newID = Identification(title,ident,z,rt,pep,prot,pepmass,premass,mods,mspre,basepeak,envelopesize,mz,intensity, homology)
                idList.append(newID)
    except:
        try: outfile.close()
        except: dummy = 0
        os.remove(out)
        return "ID"
    
    try:
        finalInts = [] # store the last intensity value of each scan to get an idea of the noise level
        count = 0 # metric of how many spectra are processed
        print 'Reading in spectrum file '+spectrumfile
        with mzml.read(spectrumfile) as reader:
            for id in reader:
                if id.get('ms level')==1:
                    count = count+1
                    mzs = list(id.get('m/z array'))
                    ints = id.get('intensity array')
                    rettime = id.get('scanList').get('scan')[0].get('scan start time')
                    try:
                        finalInts.append(ints[-1])
                    except:
                        dummy=0
    
                    # Now we get to the actual extraction part
                    # I start 1.5 minutes before and go until 1.5 minutes after
                    # the given retention time to account for drift that may
                    # occur between samples because the IDs only come from the
                    # first few samples in the time course.
                    for peptide in idList:
                        checkRT = peptide.rt
                        if abs(checkRT-rettime)<timeWindow:
                            charge = peptide.charge
                            currentIsotope = Isotope([])
                            isoArray = []
                            for i in range(0,N_ISOS):
                                peakType = 'm'+str(int(i))
                                checkmass = peptide.basepeak+(i*NEUTRON/charge)
                                reach = ppmWindow/1000000.0*checkmass
                                preIndex = bisect.bisect_left(mzs,checkmass) # binary search to find insertion point
                                # binary search doesn't always give the correct index because it tells us where to insert our peak,
                                # it doesn't give us the index of the closest peak so we have to check which slows the code down
                                minus = preIndex-1
                                plus = preIndex+1
                                try:
                                    array = [mzs[minus],mzs[preIndex],mzs[plus]]
                                except:
                                    continue
                                closestVal = min(array,key=lambda x:abs(x-checkmass))
                                try:
                                    goodIndex = mzs.index(closestVal)
                                except:
                                    goodIndex = -1
                                if abs(mzs[goodIndex]-checkmass)<reach:
                                    currentIsotope.addPeak(Peak(mzs[goodIndex],ints[goodIndex],rettime,peakType))
                                isoArray.append(ints[goodIndex])
    
                            if currentIsotope.isValid():
                                currentIsotope.setArray(isoArray[0:4])
                                peptide.addIsotope(currentIsotope)
    
        noise = 3*np.mean(finalInts)
        length = len(idList)
        failed = 0 # not enough points for good data
        passed = 0 # good data
        invalid = 0 # not found in data (should be lots of these due to extendinator.py
        print 'Finding envelopes...'
        for peptide in idList:
            isoList = peptide.getIsotope()
            goodIsos = []
            rtList = []
            angleList = []
            peakList = []
    
            # this writes out all the peaks that were previously extracted
            # the data is displayed by EnView.py
            
            minAngle = 1000000.0 # we want to subtract the smallest angle from all the rest in a sort of normalization
    
            # this loop finds the smallest angle and puts all the angles in an array
            for i in range(1,len(isoList)):
                prev = isoList[i-1].getArray()
                cur = isoList[i].getArray()
                angle = angle_between(prev,cur)
                if angle<minAngle:
                    minAngle = angle
                try:
                    angleList.append(angle)
                except: dummy = 0
    
            # now we subtract that smallest angle and antilog transform to increase spread between good and bad data
            for i in range(0,len(angleList)):
                angleList[i] = 10.0**(angleList[i]-minAngle)
    
            if len(angleList) < 3: # not enough data so we move on
                invalid = invalid + 1
                continue
    
            for i in range(0,len(angleList)): # an angle of 1-1.2 is generally good and we grab only those isotopes that correspond as such
                if angleList[i]<1.2:
                    goodIsos.append(isoList[int(i)])
    
            maxIntensity = 0 # the data gets normalized as part of the combining process so we save the maxIntensity to have an idea of the intensity of the series
            for thing in goodIsos:
                peak = thing.getM0()
                if peak.getI() > maxIntensity:
                    maxIntensity = peak.getI()
    
            for thing in goodIsos: # we must have a good enough S/N and not be tailing
                if thing.getM0().getI()>0.1*maxIntensity and thing.getM0().getI()>noise:
                    peakList.extend(thing.getPeaks())
                    rtList.append(thing.getM0().getRT())
    
            if len(rtList)<3: # not enough data so we move on
                invalid = invalid + 1
                continue
    
            # for diagnostic purposes
            minRT = -1
            maxRT = -1
            try:
                minRT = min(rtList)
                maxRT = max(rtList)
            except: dummy=0
    
            peakList = [] # now we start again and normalize  
    
            for thing in goodIsos:
                if thing.getM0().getI()>0.1*maxIntensity and thing.getM0().getI()>noise:
                    thing.normalize()
                    peakList.extend(thing.getPeaks())
    
            processedPeaks = []
    
            for i in range(0,N_ISOS):
                peakType = 'm'+str(int(i))
                try:
                    condensedPeak = condense(peakType,peakList)
                    processedPeaks.append(condensedPeak)
                except:
                    processedPeaks.append(Peak(0,0,0))
    
            
            try:
                mzvals = []
                intvals = []
                peptide.setMinRT(minRT)
                peptide.setMaxRT(maxRT)
                for entry in processedPeaks:
                    if np.isnan(entry.getMZ()):
                        mzvals.append(0.0)
                    else:
                        mzvals.append(entry.getMZ())
                    intvals.append(entry.getI()*maxIntensity)
                peptide.setSize(len(mzvals))
                peptide.setMZs(mzvals)
                peptide.setINTs(intvals)
                # if you want to pull out lots of peaks, comment out the first if and uncomment the second
                # rate calculator can't handle the nan's that will show up in unlabeled data if we go out too many peaks
                # so for now we are removing them
                #if minRT!=-1 and not (np.isnan(sum(mzvals)) or np.isinf(sum(mzvals)) or np.isnan(sum(intvals)) or np.isinf(sum(intvals))):
                if minRT!=-1:
                    outline=peptide.toString()
                    if outline != "":
                        outfile.write(outline)
                    passed = passed+1
                else: 
                    failed = failed+1
            except:
                dummy=0
    
        outfile.close()
        endTime = time.time()
        print (endTime-startTime)/60.0,' minutes.'
        print 'Processed '+str(length)+' spectra.\n'
        return ""
    except:
        try: outfile.close()
        except: dummy = 0
        os.remove(out)
        return str(spectrumfile) + ", "
