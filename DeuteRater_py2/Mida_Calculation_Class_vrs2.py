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
#basic normalization function needed for abundance
def normalize(data):
    norm_data = []
    total = sum(data)
    for x in data:
        norm_data.append(x/total)
    return norm_data
    
class mida(object):
    def __init__(self, mass, data, equations, label, unlabeled, abund, filters):
        #set the amount of labeling for future calculations
        if label == 0: self.label = filters["Zero Labeling Check"]
        else: self.label = label
        self.abund = abund
        #set the # of peaks to deal with
        if mass < filters["mass cutoff"]: size = int(filters["Peaks included if under mass cutoff"])
        else: size = int(filters["Peaks included if over mass cutoff"])
        #both abund and defect need the same data, but values are different.  poly is just the calculation that fits the coefficients, parser number is the # of coefficients
        #filter removes peaks that are not trustworthy from analysis of change (too small are untrustworthy), and true size is used for # of peaks later (defect needs an adjustment)
        if self.abund:
            self.exponential = lambda p, x: -(p[0] -p[0]*np.exp(-p[1]*x))
            self.poly = lambda p, x: p[0]*x**3 + p[1]*x**2 + p[2]*x
            self.parser_number = 3
            self.filter = filters["Minimum Abund Change"]
            self.true_size = size
        else:
            self.poly = lambda p, x: p[0]*x**4 + p[1]*x**3 + p[2]*x**2 + p[3]*x
            self.parser_number = 4
            self.true_size = size -1
        data = [float(x) for x in data] #turn the data into floats
        
        #abundance needs normalization and defect needs to be turned from m/z to distance
        if self.abund: correct_data = normalize(data[:self.true_size])
        else: 
            correct_data = []
            for r in range(1,len(data)):
                correct_data.append(data[r]-data[0])
        #extract equations
        self.equation_parser(equations[:self.true_size*self.parser_number])
        #calculate difference from baseline 
        self.delta = []
        if abund: unlabeled = normalize(unlabeled[:self.true_size])
        else: unlabeled = unlabeled[:self.true_size]
        for r in range(self.true_size):
            self.delta.append(correct_data[r]-unlabeled[r])
    #extracts the coefficients from a list of floats into coefficients for disticnt lines
    def equation_parser(self, equations):
        equations = list(equations)
        self.equations =[]
        coefficient = []
        if self.abund:
            for r in [0,1]:
                coefficient.append(equations[r])
            self.equations.append(coefficient)
            coefficient = []
            del equations[1]
            del equations[0]
        for r in range(len(equations)):
            if r%self.parser_number == 0 and r != 0:
                self.equations.append(coefficient)
                coefficient = []
            coefficient.append(equations[r])
        self.equations.append(coefficient)
    #calculate maximum change possible from the amount of labeling
    def max_peak_calc(self):
        self.max_change = []
        first = True
        for coefficients in self.equations:
            if self.abund and first:
                first = False
                self.max_change.append(self.exponential(coefficients,self.label))
            else: self.max_change.append(self.poly(coefficients,self.label))
    #function that actually calculates the change that has occured
    def change(self):
        self.max_peak_calc()
        turnover = []
        reporter = []
        if len(self.delta) >= self.true_size:
            for r in range(self.true_size):
                if self.abund:
                    if abs(self.max_change[r]) > self.filter: # if the change is too small it is error prone and so will not be used in calculations and filters.
                        turnover.append(self.delta[r]/self.max_change[r]) #actual over amount possible
                        reporter.append(self.delta[r]/self.max_change[r])
                    else:
                        reporter.append("Value {0} was not used due to low maximum possible change".format(self.delta[r]/self.max_change[r]))
                else:
                    turnover.append(self.delta[r]/self.max_change[r]) #actual over amount possible
                    reporter.append(self.delta[r]/self.max_change[r])
        return turnover, reporter