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

Store m/z, intensity, and retention time of mass spectrum peaks
"""

class Peak:
    def __init__(self,mz,i,rt,t="null"):
        self.mz = mz
        self.i = i
        self.rt = rt
        self.diff = -1
        self.t = t
        
    def getMZ(self):
        return float(self.mz)
    
    def getI(self):
        return float(self.i)
    
    def getRT(self):
        return float(self.rt)
        
    def setI(self,newI):
        self.i=newI
        
    def getDiff(self):
        return self.diff
    
    def getType(self):
        return self.t
        
    def setDiff(self,diff):
        self.diff = diff
        
    def setType(self,t):
        self.t = t

    def setString(self):
        return str(self.mz)+","+str(self.rt)+","+str(self.i)+"\n"
