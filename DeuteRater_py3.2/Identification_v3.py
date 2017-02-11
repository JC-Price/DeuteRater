# -*- coding: utf-8 -*-
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

Store information extracted form pepXML and mzML
"""

class Identification:
    def __init__(self,title,ident,charge,rt,pep,prot,pepmass,premass,mods,mspre,basepeak,envelopesize,mz,intensity, homolog):
        self.title = title
        self.ident = ident
        self.charge = float(charge)
        self.rt = float(rt)
        self.pep = pep
        self.prot = prot
        self.pepmass = pepmass
        self.premass = premass
        self.mods = mods
        self.mspre = mspre
        self.basepeak = float(basepeak)
        self.envelopesize = envelopesize
        self.mz = mz
        self.intensity = intensity
        self.minRT = -1
        self.maxRT = -1
        self.modPep = ""
        self.homolog = homolog
        changePositions = []
        qpositions = []
        if self.mods!="":
            for mod in self.mods.split(";"):
                name = mod.split("@")[0]
                if name=="Oxidation":
                    try:
                        position = int(mod.split("@")[1])-1
                        changePositions.append(position)
                    except: dummy=0
                if name=="Pyrolidone":
                    try:
                        position = int(mod.split("@")[1])-1
                        qpositions.append(position)
                    except:
                        dummy=0
                self.modPep = "".join(c.lower() if i in changePositions else c for i,c in enumerate(self.pep))
                self.modPep = "".join(c.lower() if i in qpositions else c for i,c in enumerate(self.pep))
        self.isoList = []
                
    def addIsotope(self,isotope):
        self.isoList.append(isotope)

    def getIsotope(self):
        return self.isoList
    
    def getRT(self):
        return self.rt
        
    def getMass(self):
        return self.basepeak
        
    def getPep(self):
       return self.pep
       
    def getCharge(self): 
        return self.charge

    def setRT(self,rt):
        self.rt = rt

    def setMass(self,mass):
        self.basepeak = mass

    def setMZs(self,mzs):
        self.mz = mzs

    def setINTs(self,ints):
        self.intensity = ints
        
    def setSize(self,length):
        self.envelopesize = length
    
    def setMinRT(self,val):
        self.minRT = val

    def setMaxRT(self,val):
        self.maxRT = val

    def isMatch(self,identification):
        a = identification.getPep() == self.pep
        b = identification.getCharge() == self.charge
        if a and b:
            return True
        else:
            return False
            
    def toString(self):
        mzEntry = " ".join(map(str,self.mz))
        iEntry = " ".join(map(str,self.intensity))
        returnString = '\"'+str(self.title)+'\",\"'+str(self.ident)+'\",\"'+str(self.charge)+'\",\"'+str(self.rt)+'\",\"'+str(self.pep)+'\",\"'+str(self.prot)+'\",\"'+str(self.pepmass)+'\",\"'+str(self.premass)+'\",\"'+str(self.mods)+'\",\"'+str(self.mspre)+'\",\"'+str(self.basepeak)+'\",\"'+str(self.envelopesize)+'\",\"'+str(mzEntry)+'\",\"'+str(iEntry)+'\",\"'+str(self.minRT)+'\",\"'+str(self.maxRT)+'\",\"'+str(self.homolog)+'\"\n'
        return returnString
