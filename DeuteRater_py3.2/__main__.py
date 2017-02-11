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
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import Theory_Preparation_vrs3 as tp

import Rate_Calculator_vrs8 as rc
import Grapher as g
from PyQt4 import QtCore, QtGui, uic
import Settings
import Filter_Menu
import AdvSettings
from copy import deepcopy
from ISE_class_v10a import extract

def make_folder(folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)


location =  os.path.dirname(os.path.abspath(__file__))
#$ fix default folders
Advanced_Settings = [True]
input_default= os.path.join(os.path.dirname(location), "Extracted_Files")
output_default = os.path.join(os.path.dirname(location), "Rate_Results")
user_settings = pd.Series({"N_ISOS": 5, "lowZ": 1, "hiZ": 5, "timeWindow": 1.5, "ppmWindow": 30, "Input Folder": input_default, "Output Folder": output_default,"Heavy Element": "H", "Label Key": "tissue", "Roll Up": False, "Use Abundance": True, "Use neutromer spacing": True, "Maximum Theoretical Percent": .071,"Step Size for Labeling":.001, "Minimum Non-Zero Points": 3, "Graph": "All", "asymptote": "Fixed", "fixed":1.0, "Proliferation Adjustment": 0.0000, "bias_correction": "None", "abundance_bias" : 0.0000, "neutromer spacing_bias" : 0.0000, "combined_bias" : 0.0000})
static_filters = pd.Series({"mass cutoff": 2400, "Peaks included if over mass cutoff": 5,"Peaks included if under mass cutoff": 4,  "Weight Agreement":3, "MZ Filter":0.00, "RT Filter":0.00, "Required Sequence Length" : 6, "Number of Processors": mp.cpu_count()-1, "Required Number of Labels": 10, "Minimum Abund Change":.04, "Zero Labeling Check": .05, "Abundance Agreement Filter" : .1, "neutromer spacing Agreement Filter":.1, "Combined Agreement Filter" : .1,"ZCUTOFF":1, "Error of Zero":.01, "Error of non-replicated point":.05})
#used https://www.safaribooksonline.com/blog/2014/01/22/create-basic-gui-using-pyqt/ tutorial for gui init function and form_class call
form_class = uic.loadUiType(os.path.join(location,"Main_Menu.ui"))[0]
class InteractWithUser(QtGui.QMainWindow, form_class):
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self,parent)
        self.user_settings = user_settings
        self.static_filters = static_filters
        self.setupUi(self)
        self.default_settings = deepcopy(user_settings)
        self.default_filters = deepcopy(static_filters)
        self.default_Advanced_Settings = deepcopy(Advanced_Settings)
        self.Advanced_Settings = Advanced_Settings
        self.MakeLabelingTable.clicked.connect(self.MakeTable)
        self.MIDA_Button.clicked.connect(self.MIDA_Analysis)
        self.Rate_Button.clicked.connect(self.Rate_Analysis)
        self.MichaelButton.clicked.connect(self.MichaelTest)
        self.actionExit.triggered.connect(self.ExitButton)
        self.actionChoose_Input_Folder.triggered.connect(self.choose_input_folder)
        self.actionChoose_Output_Folder.triggered.connect(self.choose_output_folder)
        self.actionChange_Settings.triggered.connect(self.Set_Settings)
        self.actionChange_Filters_2.triggered.connect(self.Set_Filters)
        self.actionAdvanced_Settings.triggered.connect(self.Adv_Settings)
        myPixmap = QtGui.QPixmap(os.path.join(location, "Logo.JPG")) # use of command from http://stackoverflow.com/questions/8687723/pyqthow-do-i-display-a-image-properly first answer accesed 5/27/2016
        self.Logo.setPixmap(myPixmap)
        #these are all to make the table interaction easier
        self.ctrl_z_primed = False
        self.ctrl_y_primed = False
        self.Labeling_Data.currentItemChanged.connect(self.store)
        self.Labeling_Data.itemChanged.connect(self.single_change) 
        QtGui.QShortcut(QtGui.QKeySequence('Ctrl+z'),self).activated.connect(self.Undo)
        QtGui.QShortcut(QtGui.QKeySequence('Ctrl+y'),self).activated.connect(self.Redo)
        QtGui.QShortcut(QtGui.QKeySequence('Ctrl+v'),self).activated.connect(self.Paste)  #http://stackoverflow.com/questions/21682261/paste-in-the-field-of-qtableview answer 1 accessed 1/28/2016
        QtGui.QShortcut(QtGui.QKeySequence('Ctrl+c'),self).activated.connect(self.Copy)
        QtGui.QShortcut(QtGui.QKeySequence('Backspace'),self).activated.connect(self.Clear_Contents)
        QtGui.QShortcut(QtGui.QKeySequence('Del'),self).activated.connect(self.Clear_Contents)
        
    #using selected indicies modeled on http://stackoverflow.com/questions/27981202/row-and-column-numbers-of-selected-cells-in-qtablewidget response 1, accessed 1/28/2016
    def Grid_Locations(self, all_data = False):
        rows = []
        columns = []
        for grid_square in self.Labeling_Data.selectedIndexes():
            rows.append(int(grid_square.row()))
            columns.append(int(grid_square.column()))
        if all_data:
            return rows, columns
        else:
            return rows[0], columns[0]
    def store(self, table_item):
        self.storage = table_item.text()
    #we're only going to allow 1 ctrl z step
    def single_change(self, table_item):
        self.Prime_Ctrl_Z([table_item.row()], [table_item.column()], True)
    def Prime_Ctrl_Z(self, rows , columns, one_item = False):
        self.ctrl_z_row = rows
        self.ctrl_z_col = columns
        if one_item:
            self.ctrl_z_data = [self.storage]
        else:
            self.ctrl_z_data = []
            for i in range(len(self.ctrl_z_row)):
                self.ctrl_z_data.append(str(self.Labeling_Data.item(self.ctrl_z_row[i], self.ctrl_z_col[i]).text()))
        self.ctrl_z_primed= True
        #after making a change, no reason for ctrl y to be active
        if self.ctrl_y_primed:
            del self.ctrl_y_row
            del self.ctrl_y_col
            del self.ctrl_y_data
            self.ctrl_y_primed = False
    def Prime_Ctrl_Y(self, rows, columns):
        #need a try except for doing before a change is made
        self.ctrl_y_row = rows
        self.ctrl_y_col = columns
        self.ctrl_y_data = []
        for i in range(len(self.ctrl_y_row)):
            self.ctrl_y_data.append(str(self.Labeling_Data.item(self.ctrl_y_row[i], self.ctrl_y_col[i]).text()))
        self.ctrl_y_primed= True
    def Undo(self):
        if self.ctrl_z_primed:
            self.Labeling_Data.blockSignals(True)
            if not self.ctrl_y_primed:
                self.Prime_Ctrl_Y(self.ctrl_z_row, self.ctrl_z_col)
            for i in range(len(self.ctrl_z_row)):
                self.Labeling_Data.item(self.ctrl_z_row[i], self.ctrl_z_col[i]).setText(self.ctrl_z_data[i])
            self.Labeling_Data.blockSignals(False)
    def Redo(self):
        if self.ctrl_y_primed:
            self.Labeling_Data.blockSignals(True)
            for i in range(len(self.ctrl_y_row)):
                self.Labeling_Data.item(self.ctrl_y_row[i], self.ctrl_y_col[i]).setText(self.ctrl_y_data[i])
            self.Labeling_Data.blockSignals(False)
    def Clear_Contents(self):
        self.Labeling_Data.blockSignals(True)
        rows, columns = self.Grid_Locations(True)
        self.Prime_Ctrl_Z(rows,columns)
        for r in range(len(rows)):
            self.Labeling_Data.item(rows[r],columns[r]).setText("")
        self.Labeling_Data.blockSignals(False)
    def Copy(self):
        #get data
        #need try except to deal with cases where the user does clicks not in the table
        try:
            rows, columns = self.Grid_Locations(True)
            min_col = min(columns) # get lowest column number
            entries = []
            col_length = max(rows)-min(rows) +1 # length (add 1 to account for the one we are on)
            for r in range(len(rows)):
                #goes by column, so first fill column list
                if columns[r] == min_col:
                    entries.append(str(self.Labeling_Data.item(rows[r], columns[r]).text()))
                #then append to that list as needed
                else:
                    current_col = columns[r]-min_col
                    entries[r - col_length*current_col] += '\t'
                    entries[r - col_length*current_col] += str(self.Labeling_Data.item(rows[r], columns[r]).text())
            #turn the list into one string that can go to the clipboard
            final  = ""
            for i in entries:
                final += i
                final +='\n'
            # from question in: http://stackoverflow.com/questions/1073550/pyqt-clipboard-doesnt-copy-to-system-clipboard accessed 1/28/2016
            clipboard = QtGui.QApplication.clipboard()
            clipboard.setText(final[:-1]) # get rid of last newline
        except:
           c = "copy"
    def Paste(self):
        #try except if they try and paste not in table.
        try:
            current_row, current_column = self.Grid_Locations() #get top left corner of user selection and start there
            text = QtGui.QApplication.clipboard().text()
            # text parsing inspired by https://riverbankcomputing.com/pipermail/pyqt/2013-May/032761.html accessed 1/28/16
            text = text.split('\n')#split rows
            text = list(text) #keep in mind we are in a pyqt gui, so all data types are Qlists and Qstrings unless you force th issue
            widest = 0
            for i in range(len(text)):
                text[i] = str(text[i]).split('\t')#split columns
                if len(list(text[i])) > widest: widest = len(text[i])
            if text[-1] == ['']: del text[-1] #excel adds whitespace to end of clipboard.  must remove 
            if len(text) <= self.Labeling_Data.rowCount() - current_row +1 and widest <= self.Labeling_Data.columnCount() - current_column:
                self.Labeling_Data.blockSignals(True)
                zrow = []
                zcol = []
                for l in range(current_row, current_row+len(text)):
                    for c in range(current_column, current_column + widest):
                        zrow.append(l)
                        zcol.append(c)
                self.Prime_Ctrl_Z(zrow, zcol)    
                for r in range(len(text)):
                    for i in range(len(text[r])):
                        self.Labeling_Data.item(current_row +r,current_column + i).setText(text[r][i])
            self.Labeling_Data.blockSignals(False)
        except:
            p = "paste"
    def get_labeling_info(self):
        #initialize needed things
        ELEMENT_LIST = ['C', 'N', 'O', 'H', 'S']
        AA_DICT = {}
        #read in the files and complain if they are not there
        try:
            aa_series = pd.read_csv(os.path.join(location, "Amino Acid Elemental Compositions.csv"))
            number_labels = pd.read_csv(os.path.join(location, "Number of Labeling Sites.csv"))
        except IOError:
            return "Amino Acid Elemental Compositions.csv or Number of Labeling Sites.csv are missing.  Correct to continue."
        #reassign the indicies
        aa_series = aa_series.set_index("Amino Acid")
        number_labels = number_labels.set_index("name")
        for amino in aa_series.index:
            if amino not in number_labels.columns:
                return "Amino acid {0} is not present in the Labeling Site Patterns.  Please correct in \"Labeling Site Patterns\" in the Advanced Settings Menu to continue.".format(amino)
            AA_DICT[amino] = []
            for elem in ELEMENT_LIST:
                AA_DICT[amino].append(aa_series[elem][amino])
            if self.user_settings["Label Key"] != "global":
                label_value = number_labels[amino][self.user_settings["Label Key"]]
                if label_value <= AA_DICT[amino][ELEMENT_LIST.index(self.user_settings["Heavy Element"])]:
                    AA_DICT[amino].append(label_value)
                else:
                    return "Labeling in labeling site pattern exceeds the number of {0} for {1}".format(self.user_settings["Heavy Element"],amino)
            else:
                AA_DICT[amino].append(AA_DICT[amino][ELEMENT_LIST.index(self.user_settings["Heavy Element"])])
        return AA_DICT
    def MichaelTest(self):
        make_folder(self.user_settings["Input Folder"])
        QtGui.QMessageBox.information(self, "Info", "Please select an id file")
        bad_files = False
        idFileName = str(QtGui.QFileDialog.getOpenFileName(self, "Select ID file", self.user_settings["Input Folder"]))
        while idFileName[-4:] != ".tsv" and idFileName[-4:] != ".csv"and idFileName != "":
             QtGui.QMessageBox.information(self, "Bad Input", "That is not a .tsv or .csv file.  Please choose a csv or tsv file.")
             idFileName = str(QtGui.QFileDialog.getOpenFileName(self, "Select ID file", self.user_settings["Input Folder"]))
        if idFileName == "":
            return
        QtGui.QMessageBox.information(self, "Info", "Please select mzml(s) to extract")   
        continue_loop = True
        if not os.path.isdir(self.user_settings["Output Folder"]):
            os.makedirs(self.user_settings["Output Folder"])
        self.user_settings.to_csv(os.path.join(self.user_settings["Output Folder"], "Settings_Extraction.csv"))
        while continue_loop:
            mzmlList = []
            bad_input = False
            data = QtGui.QFileDialog.getOpenFileNames(self, "Select mzML files", os.path.dirname(idFileName))
            if data == "": return
            else:
                for path in data:
                    if path[-5:] != ".mzML":
                        bad_input = True
                        break
                    else:mzmlList.append(str(path))
                if bad_input:
                    QtGui.QMessageBox.information(self, "Error", "One of the files you selected was not an mzML.  Please try again")
                else: continue_loop = False
        mzml_len = len(mzmlList)
        if mzml_len != 0:
            error_message = ""
            for specFile in mzmlList:
                current_pos = mzmlList.index(specFile) + 1
                self.Update_User.setText("Extracting mzML file {0} out of {1}".format(current_pos, mzml_len))
                app.processEvents()
                f = os.path.basename(specFile)
                basename = os.path.join(self.user_settings["Input Folder"], '{}.csv'.format(f.split(".")[0]))
                error = extract(idFileName,specFile,basename,self.user_settings["N_ISOS"],self.user_settings["lowZ"],self.user_settings["hiZ"],self.user_settings["timeWindow"],self.user_settings["ppmWindow"])
                if error == "ID":
                    QtGui.QMessageBox.information(self, "Error", "Your ID file has bad headers or data.  Please correct and try again.")
                    break
                elif error != "":
                    error_message += error
        try:
            if error_message != "":
                final_error = "mzML file(s) " + error_message+ "had errors in form or data and could not be extracted."
                QtGui.QMessageBox.information(self, "Error", final_error)
            elif error == "":
                QtGui.QMessageBox.information(self, "Done", "Your Extraction is Complete")
        except UnboundLocalError:
            dummy = 0
        self.Update_User.setText("Waiting for Command")
        app.processEvents()
    def choose_input_folder(self):
        make_folder(self.user_settings["Input Folder"])
        folder = QtGui.QFileDialog.getExistingDirectory(self, "Select a Folder", self.user_settings["Input Folder"], QtGui.QFileDialog.ShowDirsOnly)
        #returns an empty string if nothing is chosen so we can't have it be an empty string or stuff will error.
        if len(str(folder)) > 0:
            self.user_settings["Input Folder"] = str(folder)
    def choose_output_folder(self):
        make_folder(self.user_settings["Output Folder"])
        folder = QtGui.QFileDialog.getExistingDirectory(self, "Select a Folder",  self.user_settings["Output Folder"], QtGui.QFileDialog.ShowDirsOnly)
        if len(str(folder)) > 0:
            self.user_settings["Output Folder"] = str(folder)
    def Set_Settings(self):
        self.set_menu = Settings.Change_Settings(self, self.user_settings, self.default_settings)
        self.set_menu.show()
    def Set_Filters(self):
        self.set_menu_filter = Filter_Menu.Change_Filters(self, self.static_filters, self.default_filters)
        self.set_menu_filter.show()
    def Adv_Settings(self):
        self.set_adv_menu = AdvSettings.Change_Adv_Settings(self, self.Advanced_Settings, self.default_Advanced_Settings, self.user_settings, self.static_filters)
        self.set_adv_menu.show()    
    def ExitButton(self):
        # from http://stackoverflow.com/questions/1414781/prompt-on-exit-in-pyqt-application
        reply = QtGui.QMessageBox.question(self, "Quit Option", "Are you sure you wish to exit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            self.close()
    def MakeTable(self):
        try:
            if len(os.listdir(self.user_settings["Input Folder"])) <1 :
                #used qtgui examples
                QtGui.QMessageBox.information(self, "Empty Folder","Please choose a folder with data from the pull down menu.")
                return
        except:
            QtGui.QMessageBox.information(self, "Input Folder does not exist","Please choose a folder that exists.")
            return
        self.filenames =[]
        for thing in os.listdir(self.user_settings["Input Folder"]):
            if os.path.isfile(os.path.join(self.user_settings["Input Folder"], thing)) and thing[-4:] == ".csv":
                self.filenames.append(thing)
        if len(self.filenames) < 1:
            QtGui.QMessageBox.information(self, "Empty Folder","Folder has no .csv files.  Please choose a folder with files from the pull down menu.")
        self.Labeling_Data.blockSignals(True)
        self.Labeling_Data.setRowCount(len(self.filenames))
        for r in range(len(self.filenames)):
            self.Labeling_Data.setItem(r, 0, QtGui.QTableWidgetItem())
            self.Labeling_Data.item(r,0).setText(self.filenames[r])
            self.Labeling_Data.setItem(r, 1, QtGui.QTableWidgetItem())
            self.Labeling_Data.item(r,1).setText("")
            self.Labeling_Data.setItem(r, 2, QtGui.QTableWidgetItem())
            self.Labeling_Data.item(r,2).setText("")
        self.Labeling_Data.blockSignals(False)

    def MIDA_Analysis(self):
        #http://stackoverflow.com/questions/11956803/retrieving-data-from-columns-qtablewidget 2nd answer
        water_table = []
        too_high = False
        if not os.path.isdir(user_settings["Input Folder"]):
            QtGui.QMessageBox.information(self, "Error", "Input Folder does not exist.")
            return
        try:
            if self.Labeling_Data.rowCount() == 0:
                QtGui.QMessageBox.information(self, "Error", "Labeling Table does not exist.")
                return
            filenames = []
            for thing in os.listdir(self.user_settings["Input Folder"]):
                if os.path.isfile(os.path.join(self.user_settings["Input Folder"], thing)) and thing[-4:] == ".csv":
                   filenames.append(thing)
            if filenames != self.filenames:
                QtGui.QMessageBox.information(self, "Error", "Filenames in Input Folder do not match filenames in labeling table. Please correct before proceeding.")
                return
            for row in range(self.Labeling_Data.rowCount()):
                final_row = []
                for i in range(3):
                    thing = self.Labeling_Data.item(row, i)
                    if i == 0: 
                        if self.filenames[row] == str(thing.text()):
                            final_row.append(str(thing.text()))
                        else:
                            QtGui.QMessageBox.information(self, "Error", "Filenames have been altered.  Please correct.")
                            return
                    else: 
                        try:
                            value = float(thing.text())
                            if value < 0:
                                QtGui.QMessageBox.information(self, "Error", "Negative time or labeling make no sense.")
                                return
                            if i == 2 and value > 1:
                                QtGui.QMessageBox.information(self, "Error", "Too much label: Please enter labeling in decimal.")
                                return
                            if i == 2 and value > user_settings["Maximum Theoretical Percent"] and not Advanced_Settings[0]: too_high = True
                            final_row.append(value)
                        except:
                            QtGui.QMessageBox.information(self, "Error","Labeling Table contains text or blanks. Please correct.")
                            return
                water_table.append(final_row)
        except AttributeError:
            QtGui.QMessageBox.information(self, "Labeling Table Non-Existant","Please create table.")
        if too_high:
            message = QtGui.QMessageBox.question(self, "Quit Option", "Experimental labeling exceeds theoretical labeling.  This can cause inaccurate data.  Do you wish to continue?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if message == QtGui.QMessageBox.No:
                return
        AA_DICT = self.get_labeling_info()
        if type(AA_DICT) == str:
            QtGui.QMessageBox.information(self, "Error",AA_DICT)
            return
        print ("Beginning Analysis")
        make_folder(self.user_settings["Output Folder"])
        self.Update_User.setText("Beginning Analysis")
        app.processEvents()
        print ("Pulling in Initial Data")
        self.Update_User.setText("Pulling in Initial Data")
        app.processEvents()
        self.water = pd.DataFrame(water_table, columns = ["file", 'time', 'labeling'])
        self.water = self.water.set_index('file')
        self.water.to_csv(os.path.join(self.user_settings["Output Folder"], "Time_and_Labeling_Data.csv"))
        print ("Filtering Initial Data")
        self.Update_User.setText("Filtering Initial Data")
        app.processEvents()
        panda, errors = tp.begin_preparation(self.user_settings, self.static_filters, self.filenames, AA_DICT)
        if errors != "":
            QtGui.QMessageBox.information(self, "Error", "Files: {0} were not in the correct format or contained improper data. They will be ignored for the rest of the analysis".format(errors))
        panda = panda[panda['sequence'].str.len() >= self.static_filters["Required Sequence Length"]]
        if len(panda) == 0:
            QtGui.QMessageBox.information(self, "Insufficient IDs made it through filter",
                    "Use better data or relax the filters.")
            return
        print ("Saving Intermediate Files")
        self.Update_User.setText("Saving Intermediate Files")
        app.processEvents()
        if Advanced_Settings[0]:
            import Emass_Port_vrsFast2 as em
            import MIDA_Calculator_vrsFast2 as mc
            #fast filter
            panda = panda.join(self.water, on = "file")
            fast_panda = panda[["sequence", "labeling"]]
            fast_panda = fast_panda.drop_duplicates()
            fast_panda.to_csv(os.path.join(self.user_settings["Output Folder"], "Unique_Peptide_List.csv"), index = False, header = True)
            panda.to_csv(os.path.join(self.user_settings["Output Folder"], "Combined_Experimental_Data.csv"), index = False)
            self.user_settings.to_csv(os.path.join(self.user_settings["Output Folder"], "Settings_Isotope_Analysis.csv"))
            self.static_filters.to_csv(os.path.join(self.user_settings["Output Folder"], "Filters_Isotope_Analysis.csv"))
            fast_panda.to_csv()
            #fast_panda["labeling"] = fast_panda["labeling"].replace(to_replace = 0, value = self.static_filters["Zero Labeling Check"])
            fast_panda = fast_panda.groupby(["sequence"])
            manager = mp.Manager()
            inqueue = mp.JoinableQueue()
            for_emass = manager.dict()
            datamanip = manager.list()
            error_list = manager.list()
            maximum = 0
            for i, j in fast_panda:
                for_emass[i] = list(j["labeling"])
                maximum += len(list(j["labeling"]))
            for key in for_emass.keys():
                inqueue.put(key)
            #reset processor list
            p_list = []
            #start the processes
            num_processors = int(self.static_filters["Number of Processors"])
            if num_processors < 1:
                num_processors = 1
            for r in range(num_processors): 
                #generate function 
                p_list.append(mp.Process(target= em.emass, args = (inqueue, for_emass, datamanip, error_list, AA_DICT, self.user_settings, self.static_filters)))
                p_list[r].start() #start process
            if num_processors < mp.cpu_count():
                import time
                while not inqueue.empty():
                    a = "{0:.2f} % theory complete".format((len(datamanip)+len(error_list))/float(maximum)*100)
                    print (a)
                    self.Update_User.setText(a)
                    app.processEvents()
                    time.sleep(10)
            #wait for all the calculations to be done before proceeding
            inqueue.join()
            if len(error_list) > 0:
                error_list = pd.DataFrame(list(error_list), columns = ["Sequence", "Labeling that removed M0"])
                error_list.to_csv(os.path.join(self.user_settings["Output Folder"], "Peptides Removed for too much Label.csv"), index = False)
            if len(datamanip) > 0: 
                #Function will return a string.  This will either be a message that the Analysis succesfully completed, or that the program errored out with a description of the error  
                print ("Performing Experimental Isotope Calculations")
                self.Update_User.setText("Performing Experimental Isotope Calculations")
                app.processEvents()
                mc.Calculate(datamanip, panda, self.water, self.user_settings, self.static_filters)
                reply = QtGui.QMessageBox.question(self, "Done", "Your Analysis is Complete.  Would you like to calculate rates as well?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    self.Rate_Analysis()
                else:
                    self.Update_User.setText("Waiting for Command")
                    app.processEvents()
            else:
                QtGui.QMessageBox.information(self, "Error", "All sequences contained too much labeling to calculate")
                self.Update_User.setText("Waiting for Command")
                app.processEvents()
        else:
            import Emass_Port_vrs4 as em
            import MIDA_Calculator_vrs4 as mc
            panda.to_csv(os.path.join(self.user_settings["Output Folder"], "Combined_Experimental_Data.csv"), index = False)
            self.user_settings.to_csv(os.path.join(self.user_settings["Output Folder"], "Settings_Isotope_Analysis.csv"))
            self.static_filters.to_csv(os.path.join(self.user_settings["Output Folder"], "Filters_Isotope_Analysis.csv"))
            peptides = pd.Series(panda['sequence'].unique(), name = "sequence")
            peptides.to_csv(os.path.join(self.user_settings["Output Folder"], "Unique_Peptide_List.csv"), index = False, header = True)
            #this code runs the main theory generation through the emass program.  This portion of the code is here because windows is odd about mulitiprocessing
            #it may work in an import but this has not been tested
            #this makes the needed data stuctures.  inqueue holds the data going into multiprocessing while datamanip holds the output
            print ("Generating Theoretical Data")
            self.Update_User.setText("Generating Theoretical Data")
            app.processEvents()
            inqueue = mp.JoinableQueue()
            manager = mp.Manager()
            datamanip = manager.list()
            error_list = manager.list()
            #loads the queue with the all unique peptide sequences
            for seq in peptides:
                inqueue.put(seq)
            #reset processor list
            p_list = []
            maximum = len(peptides)
            #start the processes
            num_processors = int(self.static_filters["Number of Processors"])
            if num_processors < 1:
                num_processors = 1
            for r in range(num_processors): 
                #generate function 
                p_list.append(mp.Process(target= em.emass, args = (inqueue, datamanip, error_list, AA_DICT, self.user_settings, self.static_filters)))
                p_list[r].start() #start process
            if num_processors < mp.cpu_count():
                import time
                while not inqueue.empty():
                    a = "{0:.2f} % theory complete".format((len(datamanip)+len(error_list))/float(maximum)*100)
                    print (a)
                    self.Update_User.setText(a)
                    app.processEvents()
                    time.sleep(10)
            #wait for all the calculations to be done before proceeding
            inqueue.join()
            if len(error_list) > 0:
                error_list = pd.DataFrame(list(error_list), columns = ["Sequence", "Labeling that removed M0"])
                error_list.to_csv(os.path.join(self.user_settings["Output Folder"], "Peptides Removed for too much Label.csv"), index = False)
            if len(datamanip) > 0: 
                #Function will return a string.  This will either be a message that the Analysis succesfully completed, or that the program errored out with a description of the error  
                print ("Performing Experimental Isotope Calculations")
                self.Update_User.setText("Performing Experimental Isotope Calculations")
                app.processEvents()
                
                mc.Calculate(datamanip, panda, self.water, self.user_settings, self.static_filters)
                reply = QtGui.QMessageBox.question(self, "Done", "Your Analysis is Complete.  Would you like to calculate rates as well?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    self.Rate_Analysis()
                else:
                    self.Update_User.setText("Waiting for Command")
                    app.processEvents()
            else:
                QtGui.QMessageBox.information(self, "Error", "All sequences contained too much labeling to calculate")
                self.Update_User.setText("Waiting for Command")
                app.processEvents()
        
        
    def Rate_Analysis(self):
        if not os.path.isdir(user_settings["Output Folder"]):
            QtGui.QMessageBox.information(self, "Error", "Output Folder does not exist.")
            return 
        if os.path.isfile(os.path.join(self.user_settings["Output Folder"], "Calculation_of_Fraction_New_Protein.csv")) and os.path.isfile(os.path.join(self.user_settings["Output Folder"], "Time_and_Labeling_Data.csv")):
            columns = ["file", "Protein ID", 'Protein name', 'Experimental_Amount_of_Label (eaol) (0 = {0})'.format(self.static_filters["Zero Labeling Check"]), "sequence"]
            if user_settings["Use Abundance"]:
                columns.extend(["M0 fraction new for calculation","abundance std_dev"])
            if user_settings["Use neutromer spacing"]:
                columns.extend(["median neutromer spacing fraction new","neutromer spacing std_dev"])
            if user_settings["Use Abundance"] and user_settings["Use neutromer spacing"]:
                columns.extend(["median combined fraction new","combined std_dev"])
            try:
                complete_data = pd.read_csv(os.path.join(self.user_settings["Output Folder"], "Calculation_of_Fraction_New_Protein.csv"),  usecols = columns, low_memory = False)
            except:
                QtGui.QMessageBox.information(self, "Error", "Your Calculation_of_Fraction_New_Protein could not be read, possibly due to missing or mis-labled headers.  check you calculation type and try again.")
                return
            try:
                water = pd.read_csv(os.path.join(self.user_settings["Output Folder"], "Time_and_Labeling_Data.csv"), usecols = ["file", "time"])
                water = water.dropna()
            except:
                QtGui.QMessageBox.information(self, "Error", "Your Time_and_Labeling.csv could not be read")
                return
            if water["time"].dtype == np.float64 or water["time"].dtype == np.int64:
                for i in set(complete_data["file"]):
                    if i not in list(water["file"]):
                        QtGui.QMessageBox.information(self, "Error", "Time_and_Labeling.csv is missing file names from Calculation_of_Fraction_New_Protein or time data for them. Please correct before proceeding")
                        return 
            else:
                QtGui.QMessageBox.information(self, "Error", "Time_and_Labeling.csv has non-numbers in the time column")
                return
            water = water.set_index('file')
        else:
            QtGui.QMessageBox.information(self, "Error", "Your output is missing Calculation_of_Fraction_New_Protein or Time_and_Labeling_Data files or both.  Correct to continue.")
            return
        print ("Performing Calculations of Experimental Rates")
        self.Update_User.setText("Performing Calculations of Experimental Rates")
        app.processEvents()
        self.user_settings.to_csv(os.path.join(self.user_settings["Output Folder"], "Settings_Rate_Analysis.csv"))
        self.static_filters.to_csv(os.path.join(self.user_settings["Output Folder"], "Filters_Rate_Analysis.csv"))
        Result, Graph = rc.rate_calc(complete_data, water, self.user_settings, self.static_filters)
        for thing in Result:
            print (thing)
            if "Completed Successfully" not in thing:
                QtGui.QMessageBox.information(self, "Error", thing)
            self.Update_User.setText(thing)
            app.processEvents()
        #$ get results into an output message
        print ("Creating Graphs")
        self.Update_User.setText("Creating Graphs")
        app.processEvents()
        if self.user_settings["Graph"] != "No" and True in Graph.values():
            for analysis_type in Graph.keys():
                if self.user_settings["Graph"] == analysis_type or  self.user_settings["Graph"] == "All":
                    if Graph[analysis_type]:
                        g.Make_Graph(analysis_type, self.user_settings["Output Folder"], self.user_settings, self.static_filters)
                    elif analysis_type == "Combined" and ( not self.user_settings ["Use Abundance"] or not self.user_settings ["Use neutromer spacing"]): 
                        print ("{0} graphs were not possible due to user settings.".format(analysis_type))
                    elif (analysis_type == "Abundance" and not self.user_settings ["Use Abundance"]) or (analysis_type == "neutromer spacing" and not self.user_settings ["Use neutromer spacing"]):
                        print ("{0} graphs was not possible due to user settings.".format(analysis_type))
                    else:
                        QtGui.QMessageBox.information(self, "Error",
                        "Graphs for {0} could not be completed due to lack of data".format(analysis_type))
        QtGui.QMessageBox.information(self, "Done", "Your Analysis is Complete")
        self.Update_User.setText("Waiting for Command")
        app.processEvents()
        
    def Magic_Button(self):
        import MIDA_Calculator_vrs4 as mc
        message = mc.Recalculate(self.user_settings, self.static_filters)
        QtGui.QMessageBox.information(self, "Done", message)
                        
if __name__ == '__main__':
    import sys
    mp.freeze_support()
    app = QtGui.QApplication(sys.argv)
    cheese = InteractWithUser(None)
    cheese.show()
    app.exec_()
