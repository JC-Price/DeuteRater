# -*- coding: utf-8 -*-
"""
Copyright (c) 2016 Michael, Porter, Bradley Naylor, J.C. Price, and Brigham Young University

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
import Peak
from Peak import *

class Isotope:
    def __init__(self,peakList):
        self.peakList = peakList
        self.defectTheta = 0
        self.abundTheta = 0
        self.array = []
        
    def addPeak(self,peak):
        self.peakList.append(peak)
        
    def removePeak(self,peak):
        self.peakList.remove(peak)
        
    def normalize(self):
        #self.peakList.sort(key=lambda x: x.getMZ())
        factor = 1
        for thing in self.peakList:
            if thing.getType()=='m0':
                factor = thing.getI()
        for thing in self.peakList:
            newI = thing.getI()/factor
            thing.setI(newI)
    
    def getPeaks(self):
        return self.peakList
        
    def getSize(self):
        return len(self.peakList)

    def hasM0(self):
        for thing in self.peakList:
            if thing.getType()=='m0':
                return True
        return False

    def getM0(self):
        for thing in self.peakList:
            if thing.getType()=='m0':
                return thing
        return None

    def isValid(self):
        a=False
        b=False
        c=False
        d=False
        for thing in self.peakList:
            if thing.getType()=='m0':
                a = True
            if thing.getType()=='m1':
                b = True
            if thing.getType()=='m2':
                c = True
            if thing.getType()=='m3':
                d = True
        return (a and b and c and d)

    def setDefectAngle(self,angle):
        self.defectTheta = angle

    def getDefectAngle(self):
        return self.defectTheta

    def setAbundAngle(self,angle):
        self.abundTheta = angle

    def getAbundAngle(self):
        return self.abundTheta

    def setArray(self,array):
        self.array = array

    def getArray(self):
        return self.array
