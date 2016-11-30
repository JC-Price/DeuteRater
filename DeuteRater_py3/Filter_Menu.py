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
from PyQt4 import QtCore, QtGui, uic
import os
import multiprocessing as mp
location =  os.path.dirname(os.path.abspath(__file__))
form_class = uic.loadUiType(os.path.join(location,"Filters_Menu.ui"))[0]
#set max for processors when this opens
# when assinging it a value check that it is between 1 and max inclusive
class Change_Filters(QtGui.QDialog, form_class):
    def __init__(self, parent = None, stuff = None, defaults = None):
        super(Change_Filters, self).__init__(parent)
        self.filters = stuff #make the settings passed in a self variable for easier manipulation (this will pass back to the main without a return)
        self.setWindowTitle("Filter Menu")
        self.setupUi(self)
        self.Processors.setMaximum(mp.cpu_count())
        self.defaults = self.make_adjustments(defaults)# get the default values for the default reset button.
        self.last_save = self.make_adjustments(stuff)
        self.adjust_values(self.last_save)
        self.ExitMenu.clicked.connect(self.close)
        self.SaveFilters.clicked.connect(self.getdata)
        self.ResetToDefaults.clicked.connect(self.default_button)
        self.Last_Saved_Filters.clicked.connect(self.saved_button)
    def default_button(self):
        self.adjust_values(self.defaults)
    def saved_button(self):
        self.adjust_values(self.last_save)
    def adjust_values(self, dictionary):
        for key in dictionary.keys():
            key.setValue(dictionary[key])
    #retrieve default settings
    def make_adjustments (self, thing):
        big_dict = {self.masscutoff:thing["mass cutoff"], self.overmasscutoff:thing["Peaks included if over mass cutoff"], self.undermasscutoff:thing["Peaks included if under mass cutoff"],  self.theory_mass_diff:thing["Weight Agreement"], self.mzfilter:thing["MZ Filter"], self.RTfilter:thing["RT Filter"], self.SeqLength:thing["Required Sequence Length"], self.Processors:thing["Number of Processors"], self.RequiredLabels:thing["Required Number of Labels"], self.minabundchange:thing["Minimum Abund Change"], self.zero_label:thing["Zero Labeling Check"], self.Abundstddev:thing["Abundance Agreement Filter"], self.neutromerstddev:thing["neutromer spacing Agreement Filter"], self.combindedstddev:thing["Combined Agreement Filter"],self.Zcutoff:thing["ZCUTOFF"], self.errorofzero:thing["Error of Zero"], self.nonreppoint:thing["Error of non-replicated point"]}
        if big_dict[self.zero_label] < 1:
            big_dict[self.zero_label] = int(big_dict[self.zero_label]*100) # get into percentage
        if big_dict[self.Processors] <1:
            big_dict[self.Processors] = 1 # in the main there is a possibility for zero and we don't want errors
        return big_dict
    def getdata(self):
        if int(self.overmasscutoff.value()) < int(self.undermasscutoff.value()):
            QtGui.QMessageBox.information(self, "Error", "Peaks if under mass cutoff cannot exceed peaks if over mass cutoff")
            return
        if int(self.Processors.value()) == mp.cpu_count() and mp.cpu_count != 1:
            reply = QtGui.QMessageBox.question(self, "All Processors", "You have selected to use all processors.  This will speed up theory caclulations, but slow down any other processes on the computer and disable percentage updates during theoretical calculation.  Is this correct?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return
        self.filters["mass cutoff"] = int(self.masscutoff.value())
        self.filters["Peaks included if over mass cutoff"] = int(self.overmasscutoff.value())
        self.filters["Peaks included if under mass cutoff"] = int(self.undermasscutoff.value())
        self.filters["Weight Agreement"] = float(self.theory_mass_diff.value())
        self.filters["MZ Filter"] = float(self.mzfilter.value())
        self.filters["RT Filter"] = float(self.RTfilter.value())
        self.filters["Required Sequence Length"] = int(self.SeqLength.value())
        self.filters["Number of Processors"] = int(self.Processors.value())
        self.filters["Required Number of Labels"] = int(self.RequiredLabels.value())
        self.filters["Minimum Abund Change"] = float(self.minabundchange.value())
        self.filters["Zero Labeling Check"] = float(self.zero_label.value())/100.0
        self.filters["Abundance Agreement Filter"] = float(self.Abundstddev.value())
        self.filters["neutromer spacing Agreement Filter"] = float(self.neutromerstddev.value())
        self.filters["Combined Agreement Filter"] = float(self.combindedstddev.value())
        self.filters["ZCUTOFF"] = float(self.Zcutoff.value())
        self.filters["Error of Zero"] = float(self.errorofzero.value())
        self.filters["Error of non-replicated point"] = float(self.nonreppoint.value())
        self.last_save = self.make_adjustments(self.filters)