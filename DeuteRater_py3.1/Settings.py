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
import sys,os
import csv
location =  os.path.dirname(os.path.abspath(__file__))
form_class = uic.loadUiType(os.path.join(location,"Settings_Menu.ui"))[0]

#make the labeling_conditions menu
labeling_conditions = []
with open (os.path.join(location,"Number of Labeling Sites.csv"), "rt") as infile:
    reader = csv.reader(infile)
    next(reader, None)
    for row in reader:
        labeling_conditions.append("{0} ({1})".format(row[0], row[1]))
labeling_conditions.append('global') #add special values
class Change_Settings(QtGui.QDialog, form_class):
    def __init__(self, parent = None, stuff = None, defaults = None):
        super(Change_Settings, self).__init__(parent)
        self.settings = stuff #make the settings passed in a self variable for easier manipulation (this will pass back to the main without a return)
        self.setWindowTitle("Settings Menu")
        self.setupUi(self)
        """deal with label key"""
        #set the labeling pattern combo box
        #attach the buttons to their functions.
        self.LabelPattern.clear()
        self.LabelPattern.addItems(labeling_conditions)
        self.defaults = self.make_adjustments(defaults)# get the default values for the default reset button.
        self.last_save = self.make_adjustments(stuff)
        self.adjust_values(self.last_save)
        self.ExitMenu.clicked.connect(self.close)
        self.SaveSettings.clicked.connect(self.getdata)
        self.ResetToDefaults.clicked.connect(self.default_button)
        self.Last_Saved_Settings.clicked.connect(self.saved_button)
    def default_button(self):
        self.adjust_values(self.defaults)
    def saved_button(self):
        self.adjust_values(self.last_save)
    #sets all values to their default settings
    def adjust_values(self, dictionary):
        for key in dictionary.keys(): 
            try:
                index = key.findText(dictionary[key])
                key.setCurrentIndex(index)
            except:
                key.setValue(dictionary[key])
    #retrieve default settings
    def make_adjustments (self, thing):
        big_dict = {self.n_isos:thing["N_ISOS"], self.lowz:thing["lowZ"], self.highz:thing["hiZ"], self.time_window:thing["timeWindow"], self.ppm_window:thing["ppmWindow"], self.HeavyElement:thing["Heavy Element"],self.LabelPattern:thing["Label Key"]}
        if thing["Roll Up"]:
            big_dict[self.RollUp] = "Yes"
        else: big_dict[self.RollUp] = "No"
        if thing["Use Abundance"] and thing["Use neutromer spacing"]:big_dict[self.CalculationMethodOptions] = "Both"
        elif thing["Use Abundance"]: big_dict[self.CalculationMethodOptions] = "Abundance"
        else: big_dict[self.CalculationMethodOptions] = "Neutromer Spacing"
        big_dict[self.TheoryLabel] =thing["Maximum Theoretical Percent"] *100
        big_dict[self.StepSize] = thing["Step Size for Labeling"] *100
        big_dict[self.TimePoints] = thing["Minimum Non-Zero Points"]
        big_dict[self.GraphOptions] = thing["Graph"]
        big_dict[self.asymptote_calc] = thing["asymptote"]
        big_dict[self.FAValue] = thing["fixed"]
        big_dict[self.Prolif_Adjust] = thing["Proliferation Adjustment"]
        big_dict[self.bias_selection]= thing["bias_correction"]
        big_dict[self.intensity_bias_box] = thing["abundance_bias"]
        big_dict[self.spacing_bias_box] = thing["neutromer spacing_bias"]
        big_dict[self.combined_bias_box] = thing["combined_bias"]
        if thing["Label Key"] != "global":
            big_dict[self.LabelPattern] = "{0} ({1})".format(big_dict[self.LabelPattern], big_dict[self.HeavyElement])
        return big_dict
    #retrieve data from user
    def getdata(self):
        if int(self.lowz.text()) > int(self.highz.text()):
            QtGui.QMessageBox.information(self, "Error", "Low charge cannot be greater than high charge!")
            return
        if str(self.LabelPattern.currentText()) != "global" and str(self.LabelPattern.currentText())[-2] != str(self.HeavyElement.currentText()):
            QtGui.QMessageBox.information(self, "Error", "Labeling Sites must match chosen Heavy Element")
            return
        if str(self.LabelPattern.currentText()) != "global":self.settings["Label Key"] = str(self.LabelPattern.currentText())[:-4]
        else:self.settings["Label Key"] = str(self.LabelPattern.currentText())
        #these ones can't be given bad values from the user interface so we don't need to check them
        self.settings["Heavy Element"] = str(self.HeavyElement.currentText())
        if str(self.RollUp.currentText())== "Yes":
            self.settings["Roll Up"] = True
        else:
            self.settings["Roll Up"] = False
        if str(self.CalculationMethodOptions.currentText()) == "Both" or str(self.CalculationMethodOptions.currentText()) == "Abundance":
            self.settings["Use Abundance"] = True
        else:
            self.settings["Use Abundance"] = False
        if str(self.CalculationMethodOptions.currentText()) == "Both" or str(self.CalculationMethodOptions.currentText()) == "Neutromer Spacing":
            self.settings["Use neutromer spacing"] = True
        else:
            self.settings["Use neutromer spacing"] = False
        self.settings["Maximum Theoretical Percent"] = float(self.TheoryLabel.value())/100+.001
        self.settings["Step Size for Labeling"] = float(self.StepSize.value())/100
        self.settings["Minimum Non-Zero Points"] = int(self.TimePoints.value())
        self.settings["Graph"] = str(self.GraphOptions.currentText())
        self.settings["N_ISOS"] = int(self.n_isos.value())
        self.settings["lowZ"] = int(self.lowz.value())
        self.settings["hiZ"] = int(self.highz.value())
        self.settings["timeWindow"] = float(self.time_window.value())
        self.settings["ppmWindow"] = float(self.ppm_window.value())
        self.settings["asymptote"] =str(self.asymptote_calc.currentText())
        self.settings["fixed"] = float(self.FAValue.value())
        self.settings["Proliferation Adjustment"] = float(self.Prolif_Adjust.value())
        self.settings["bias_correction"] = str(self.bias_selection.currentText())
        self.settings["abundance_bias"]=float(self.intensity_bias_box.value())
        self.settings["neutromer spacing_bias"]=float(self.spacing_bias_box.value())
        self.settings["combined_bias"]=float(self.combined_bias_box.value())
        self.last_save = self.make_adjustments(self.settings)
#this is for running code independantly (for testing purposes)
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = Change_Settings()
    myapp.show()
    sys.exit(app.exec_())