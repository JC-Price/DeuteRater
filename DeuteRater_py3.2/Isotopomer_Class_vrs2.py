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
ELEMENT_LIST = ['C', 'N', 'O', 'H', 'S']
class peptide:
    def __init__(self, sequence, label, AA_DICT):
        self.sequence = sequence
        self.label = label
        self.AA_DICT = AA_DICT
        self.elem_totals = {'C':0, 'N':0, 'O':1, 'H':2, 'S':0, 'X':0} #tail ends of the peptide
        self.elemental_comp()
    def Composition(self, emass=False):
        composition = ''
        for element in ELEMENT_LIST:
            if element == self.label and emass:
                composition += (element + str(self.elem_totals[element] - self.get_n()))
                composition += ('X' + str(self.get_n()))
            elif element != 'X':
                composition += (element + str(self.elem_totals[element]))
        return composition
    def MW (self):
        #standard atomic weights from http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
        #this only works for no heavy isotopes
        mw = self.elem_totals['C'] * 12.0
        mw += self.elem_totals['H'] * 1.00782503207
        mw += self.elem_totals['N'] * 14.0030740048
        mw += self.elem_totals['O'] * 15.99491461956
        mw += self.elem_totals['S'] * 31.972071 
        return mw
    def get_n (self):
        return int(self.elem_totals['X'] + .5)
    def elemental_comp(self):
        for AA in self.sequence:
            self.elem_totals['C'] += self.AA_DICT[AA][ELEMENT_LIST.index("C")]
            self.elem_totals['N'] += self.AA_DICT[AA][ELEMENT_LIST.index("N")]
            self.elem_totals['O'] += self.AA_DICT[AA][ELEMENT_LIST.index("O")]
            self.elem_totals['H'] += self.AA_DICT[AA][ELEMENT_LIST.index("H")]
            self.elem_totals['S'] += self.AA_DICT[AA][ELEMENT_LIST.index("S")]
            self.elem_totals['X'] += self.AA_DICT[AA][-1]