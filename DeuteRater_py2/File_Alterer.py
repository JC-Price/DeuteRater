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

#am not currently familiar enough with inheritance using pyqt to mess with the inheritance
#this could seriously benefit from inheritance so change that as time permits

from PyQt4 import QtCore, QtGui, uic
import sys,os
import csv
import pandas as pd

location =  os.path.dirname(os.path.abspath(__file__))

elem_comp_class = uic.loadUiType(os.path.join(location,"Elem_Comp.ui"))[0]
poss_label_class = uic.loadUiType(os.path.join(location,"Labeling_Sites.ui"))[0]
mod_class = uic.loadUiType(os.path.join(location,"Modifications.ui"))[0]

#$ deletion mechanic for elem comp and poss label

def reset_table(table, filename):
    table.setRowCount(0) #should clear? http://stackoverflow.com/questions/29266452/how-can-i-clear-a-pyqt-qtablewidget
    make_table(table, filename)

def make_table(table, filename):
    data = []
    with open(filename, 'rb') as infile:
        reader = csv.reader(infile)
        header = next(reader)
        for row in reader:
            data.append(row)
    table.setRowCount(len(data))
    table.setColumnCount(len(data[0]))
    table.setHorizontalHeaderLabels(header)
    for d in range(len(data)):
        for i in range(len(data[d])):
           table.setItem(d, i, QtGui.QTableWidgetItem())
           table.item(d,i).setText(data[d][i])    

def additional_row(table):
    new_count = table.rowCount() + 1
    table.setRowCount(new_count) #make sure it doesn't clear anything
    for i in range(table.columnCount()):
        table.setItem(new_count-1, i, QtGui.QTableWidgetItem())
        table.item(new_count-1,i).setText("")  
        
class alter_element_composition(QtGui.QDialog, elem_comp_class):
    def __init__(self, parent = None):
        super(alter_element_composition, self).__init__(parent)
        self.setupUi(self)  
        self.setWindowTitle("Elemental Composition Menu")
        self.filename =  os.path.join(location,"Amino Acid Elemental Compositions.csv")
        make_table(self.ElemTable,self.filename)
        self.save.clicked.connect(self.save_table)
        self.reset.clicked.connect(lambda: reset_table(self.ElemTable, self.filename))
        self.new_row.clicked.connect(lambda: additional_row(self.ElemTable))
        self.exit.clicked.connect(self.close)
        
    def save_table(self):
        header = []
        codes = []
        for i in range(self.ElemTable.columnCount()):
            header.append(self.ElemTable.horizontalHeaderItem(i).text())
        final = [header]
        for r in range(self.ElemTable.rowCount()):
            row = []
            delete_row = False
            code = self.ElemTable.item(r,0).text()
            description = self.ElemTable.item(r,1).text() 
            if code == "" and description =="":
                delete_row = True
            elif code == "" or description =="":
                QtGui.QMessageBox.information(self, "Error", "Amino Acid or description is blank in row {0}.  fill blanks or delete all cells in the row to save".format(r+1)) #functional
                return
            if len(code) > 1:
                QtGui.QMessageBox.information(self, "Error", "Amino Acid can only be one character in length. {0} is too long. Correct to save".format(code)) #functional
                return
            if code in codes:
                QtGui.QMessageBox.information(self, "Error", "{0} is in Amino Acid column twice.  This is not allowed. Correct to save".format(code))#functional
                return
            codes.append(code)
            row.extend([code, description])
            for i in range(2, len(header)):
                try:
                    number = str(self.ElemTable.item(r,i).text()) 
                    if delete_row and number != "":
                        QtGui.QMessageBox.information(self, "Error", "blanks in the table.  fill blanks or delete all cells in the row to save") # functional
                        return
                    elif not delete_row:
                        if not float(number).is_integer():
                            QtGui.QMessageBox.information(self, "Error", "One of the values in column {0} is not an integer.  correct to save".format(header[i])) #functional
                            return
                        if int(number) <0:
                            QtGui.QMessageBox.information(self, "Error", "One of the values in column {0} is negative.  correct to save".format(header[i])) #functional
                            return
                        row.append(int(number))
                except ValueError:
                    QtGui.QMessageBox.information(self, "Error", "One of the values in column {0} cannot be turned into a number.  correct to save".format(header[i]))#functional
                    return
            if not delete_row:
                final.append(row)
        with open(self.filename, 'wb') as outfile:
            writer = csv.writer(outfile)
            for row in final:
                writer.writerow(row)
        reset_table(self.ElemTable, self.filename)
        
class alter_mods(QtGui.QDialog, mod_class):
    def __init__(self, parent = None):
        super(alter_mods, self).__init__(parent)
        self.filename = os.path.join(location,"Modifications.csv")
        self.setupUi(self)
        self.setWindowTitle("Modifications Menu")
        make_table(self.modTable,self.filename)
        self.save.clicked.connect(self.save_table)
        self.reset.clicked.connect(lambda: reset_table(self.modTable, self.filename))
        self.add_row.clicked.connect(lambda: additional_row(self.modTable))
        self.exit.clicked.connect(self.close)
    def save_table(self): # check how this responds to header
        header = []
        for i in range(self.modTable.columnCount()):
            header.append(self.modTable.horizontalHeaderItem(i).text())
        final = [header]
        codes = []
        mods = []
        for r in range(self.modTable.rowCount()): 
            mod = self.modTable.item(r,0).text()
            code = self.modTable.item(r,1).text()
            if (code == "" and mod != "") or (code != "" and mod == ""):  
                QtGui.QMessageBox.information(self, "Error", "A row contains blanks.  please delete all cells or fill in the blanks in order to save")#functional
                return
            if code != "" and len(code) != 1:
                QtGui.QMessageBox.information(self, "Error", "Letter code must be only a single character in length.  correct to save")#functional
                return
            if " " in mod or "\t" in mod:
                QtGui.QMessageBox.information(self, "Error", "no spaces or tabs in modification names.  correct to save")#functional
                return
            if code in codes:
                QtGui.QMessageBox.information(self, "Error", "You have replicate one letter codes.  This is not allowed.  correct to save") #functional
                return
            if mod in mods:
                QtGui.QMessageBox.information(self, "Error", "You have replicate Modification names.  This is not allowed.  correct to save")#functional
                return
            if code != "" and mod != "":
                final.append([mod, code])
                codes.append(code)
                mods.append(mod)
        with open(self.filename, 'wb') as outfile:
            writer = csv.writer(outfile)
            for row in final:
                writer.writerow(row)
        reset_table(self.modTable, self.filename)
            
class amount_of_labels(QtGui.QDialog, poss_label_class):
    def __init__(self, parent = None):
        super(amount_of_labels, self).__init__(parent)
        self.filename = os.path.join(location,"Number of Labeling Sites.csv")
        data = pd.read_csv(self.filename)
        self.data = data.set_index("name")
        for i in self.data.columns:
            if i != "label":
                self.data[i] = self.data[i].astype(float)
        self.check = pd.read_csv(os.path.join(location,"Amino Acid Elemental Compositions.csv"))
        self.AAs = list(self.data.columns[1:])  #$ need some form of check that AAs are present in each.  will come back to this.
        self.setupUi(self)
        self.setWindowTitle("Labeling Menu")
        self.elems = [self.heavy_element.itemText(i) for i in range(self.heavy_element.count())] # http://stackoverflow.com/questions/7479915/getting-all-items-of-qcombobox-pyqt4-python
        self.deal_with_AAs()
        self.set_pattern_pulldown()
        self.exit.clicked.connect(self.close)
        self.Load.clicked.connect(self.load_pattern)
        self.save.clicked.connect(self.save_pattern)
        self.load_pattern()
    def set_pattern_pulldown(self, current_name = ""):
        patterns = [""]
        patterns.extend(self.data.index)
        self.PatternName.clear()
        self.PatternName.addItems(patterns)
        if current_name != "":
            i = patterns.index(current_name)
            self.PatternName.setCurrentIndex(i)
        
    def load_pattern(self):
        pattern = str(self.PatternName.currentText())
        self.Save_Name.setText(pattern)
        if pattern == "":
            self.heavy_element.setCurrentIndex(0)
        else:
            label = self.data["label"][pattern]
            self.heavy_element.setCurrentIndex(self.elems.index(label))
        self.LabelTable.setRowCount(0) #clear the table completely
        self.LabelTable.setRowCount(len(self.AAs))
        self.LabelTable.setColumnCount(2)
        self.LabelTable.setHorizontalHeaderLabels(["Amino Acid", "Possible Labeling Sites"])
        for r in range(len(self.AAs)):
            self.LabelTable.setItem(r,0, QtGui.QTableWidgetItem())
            self.LabelTable.setItem(r,1, QtGui.QTableWidgetItem())
            self.LabelTable.item(r,0).setText(self.AAs[r])
            if pattern == "": self.LabelTable.item(r,1).setText("0")
            else: self.LabelTable.item(r,1).setText(str(self.data[self.AAs[r]][pattern]))
        
    def save_pattern(self):
        name = str(self.Save_Name.text())
        completed_AAs = []
        if name == "":
            QtGui.QMessageBox.information(self, "Error", "New Name must be present to save.  Correct and try again.")
            return
        replace =False
        if name in self.data.index:
            reply = QtGui.QMessageBox.question(self, "Name exists", "Save name already exists.  Do you really want to overwrite it?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                QtGui.QMessageBox.information(self, "Info", "Data was not saved")
                return
            else: replace = True
        new_data = pd.Series()
        new_data.name = name
        new_data["label"] = str(self.heavy_element.currentText())
        for r in range(self.LabelTable.rowCount()):
            amino_acid_code = str(self.LabelTable.item(r,0).text())
            X_number = self.LabelTable.item(r,1).text()
            #no blanks (other tables can delete not here)
            #functional
            if amino_acid_code == "" or X_number == "":
                QtGui.QMessageBox.information(self, "Error", "Table contains blanks. please correct to save")
                return
            #needs to be a valid AA (something we have an elemental composition for)
            #functional
            if amino_acid_code not in self.AAs:
                QtGui.QMessageBox.information(self, "Error", "amino acid {} is not a valid Amino Acid. please correct to save.".format(amino_acid_code))
                return
            #no repeats
            #functional
            if amino_acid_code in completed_AAs:
                QtGui.QMessageBox.information(self, "Error", "amino acid {} is present in the table twice. please correct to save.".format(amino_acid_code))
                return
            else:
                completed_AAs.append(amino_acid_code)
            #can only be one character long, but no need to check that because that should be taken care of in AA page
            #test if number is a numerical value between 0 and max number for that element
            try:
                #functional
                X_number = float(X_number)
                #these next two lines are inefficient.  consider revising.
                check_value = self.check[self.check["Amino Acid"]==amino_acid_code]
                upper_limit = check_value[new_data["label"]].values[0]
                if X_number <0 or X_number >upper_limit:
                    QtGui.QMessageBox.information(self, "Error", "amino acid {0} has labeling that is not between 0 and the amount of {1} for that amino acid ({2}). please correct to save.".format(amino_acid_code,self.heavy_element.currentText(),upper_limit))
                    return
            except ValueError:
                #functional
                QtGui.QMessageBox.information(self, "Error", "amino acid {}\'s potential labels is not a number. please correct to save.".format(amino_acid_code))
                return
            new_data[amino_acid_code] = X_number
        if replace:
            for aa_code in new_data.iteritems():
                self.data.set_value(name, aa_code[0], aa_code[1])
        else:
            self.data = self.data.append(new_data)
            self.set_pattern_pulldown(name)
        self.data.to_csv(self.filename)
    def deal_with_AAs(self):
        #functions
        for a in self.check["Amino Acid"]:
            if a not in self.AAs:
                zeroes = []
                for i in range(len(self.data.index)):
                    zeroes.append(0.0)
                self.data[a] = zeroes
                self.AAs.append(a)
        #functions
        delete = []
        for a in self.AAs:
            if a not in self.check["Amino Acid"].values:
                 delete.append(a)
        if len(delete) > 0:
            reply = QtGui.QMessageBox.question(self, "Extra Labeling", "The Following Amino Acids have potential labels, but no elemental compositions: {0}.  Would you like to delete them from labeling?".format(delete), QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                for d in delete:
                    del self.data[d]
                    i = self.AAs.index(d)
                    del(self.AAs[i])
            
            