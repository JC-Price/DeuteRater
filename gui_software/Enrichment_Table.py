# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Christian Andersen, Chad Quilling, J.C. Price, and Brigham Young University
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
"""
this governs the table the user fills in to provide the enrichment for a given subject

"""

import os
import pandas as pd

from PyQt5 import uic, QtWidgets, QtCore, QtGui

import gui_software.general_table_functions as gtf


#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Enrichment_Table.ui")

loaded_ui = uic.loadUiType(ui_file)[0]

#  it is tricky to actually get the header out of the qtablewidget and they
#  need different checks anyway so we'll just declare it here
current_columns = ["Subject ID", "Time", "Enrichment"]

class EnrichmentWindow (QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, min_allowed_times = None, 
                 starting_num_timepoints = None, outfile = None,
                 max_enrichment = None):
        super(EnrichmentWindow, self).__init__(parent)
        # allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        self.early_exit = True
        self.setupUi(self)
        
        self.max_enrichment = max_enrichment
        self.outfile = outfile
        self.min_allowed_times = min_allowed_times
        temp_df = pd.read_csv(outfile, sep = "\t")
        self.subject_ids = list(temp_df[current_columns[0]].unique())
        gtf.set_up_table(self.EnrichmentTable, self.subject_ids, 
                         starting_num_timepoints)
        
        #  just set the minimum from a setting to avoid needing multiple
        # if statements
        self.enrichments_for_subject.setMinimum(self.min_allowed_times)
        self.enrichments_for_subject.setValue(starting_num_timepoints)
        self.current_timepoints = starting_num_timepoints
        
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)
        self.UpdateTableButton.clicked.connect(self.update_table_size)
        self.setWindowTitle("Fill in Enrichmet Values")
        
        # make some easy table keyboard shortcuts. 
        # only do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(lambda: gtf.Copy_to_Clipboard(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(lambda:gtf.Paste(
                                self.EnrichmentTable))

    # allows increasing or decreasing table size
    def update_table_size(self):
        desired_value = int(self.enrichments_for_subject.value())
        # no point in doing anything if not adding or deleting
        if desired_value == self.current_timepoints:
            return
        elif desired_value < self.current_timepoints:
            row_to_delete = []
            for s in range(len(self.subject_ids)):
                #  offset finds the indicies of the first instance of the 
                # subject id
                offset = self.current_timepoints * s
                row_to_delete.extend([offset + x for x in range(desired_value,
                                              self.current_timepoints)])
            # go backwards to ensure that we don't move critical rows into 
            # delete positions
            for index_value in reversed(row_to_delete):
                self.EnrichmentTable.removeRow(index_value)
        else:
            for s in range(len(self.subject_ids)):
                offset = desired_value * s
                for index_value in range(self.current_timepoints, desired_value):
                    true_index = index_value + offset
                    self.EnrichmentTable.insertRow(true_index)
                    name_item = QtWidgets.QTableWidgetItem(self.subject_ids[s])
                    name_item.setFlags(QtCore.Qt.ItemIsSelectable | 
                                   QtCore.Qt.ItemIsEnabled)
                    self.EnrichmentTable.setItem(true_index, 0, name_item)
                    for i in range(1, self.EnrichmentTable.columnCount()):
                        # makes an empty string in all other cells
                        self.EnrichmentTable.setItem(true_index, i, 
                                        QtWidgets.QTableWidgetItem(""))
        self.current_timepoints = desired_value
    
    # check that we may save and continue.  if not complain to the user and give
    # don't continue so they can correct the problem
    def check_and_save (self):
        # don't have to check the ids, those are out of the user's control in this table
        class_dict = {}
        for r in range(self.EnrichmentTable.rowCount()):
            # first column is not editable so can be ignored
            current_id = str(self.EnrichmentTable.item(r,0).text())
            #  can't just loop as with Time_Table.py since we need to analyze
            #  together and need to store in a class
            time_string = self.EnrichmentTable.item(r,1).text()
            enrichment_string = self.EnrichmentTable.item(r,2).text()
            # blanks are okay as long as all are blank
            if  time_string == "" and enrichment_string == "":
                continue
            time_value = gtf.basic_number_error_check(
                    time_string, current_columns[1], r)
            # need to block enrichment over 100%
            enrichment_value =  gtf.basic_number_error_check(
                    enrichment_string, current_columns[2], r, self.max_enrichment, True)
            for value in [time_value, enrichment_value]: 
                if type(value) == str:
                    QtWidgets.QMessageBox.information(self, "Error", value)
                    return 
            if current_id in class_dict.keys():
                class_dict[current_id].add_data(time_value, enrichment_value)
            else:
                class_dict[current_id] = subject(current_id, time_value,
                                                    enrichment_value)
        # if all are blank (or all for one subject are blank) it will hit
        # continue and never be added.  this is somehting to complain about
        for id_value in self.subject_ids:
            if id_value not in class_dict.keys():
                QtWidgets.QMessageBox.information(self, "Error", ("Subject ID {}"
                            " has no values. Correct to proceed".format(id_value)))
                return 
        
        # now that we have the data in the subject object and let it do the 
        # error checking
        for sample_id in class_dict.keys():
            return_value = class_dict[sample_id].error_checking(
                self.min_allowed_times)
            if type(return_value) == str:
                QtWidgets.QMessageBox.information(self, "Error", return_value)
                return
        self.make_save_file(class_dict)
            
    def make_save_file(self, class_dict):
        # csv would be faster but since the files may be of different lengths
        # it would likely be more complex than necessary
        full_names = []
        full_times = []
        full_enrichment = []
        blank_column = [""]
        for key in class_dict.keys():
            full_names.extend([key for i in range(len(class_dict[key].times))])
            full_times.extend(class_dict[key].times)
            full_enrichment.extend(class_dict[key].enrichments)
        blank_column = blank_column * len(full_names)
        df = pd.read_csv(self.outfile, delimiter ="\t")
        df2 = pd.DataFrame({"Spacing Column 1": blank_column, 
                            "Subject ID Enrichment":full_names,
                            "Time Enrichment": full_times,
                            "Enrichment": full_enrichment
                            })
        df = pd.concat([df, df2], axis =1)
        df.to_csv(self.outfile, sep ="\t", index = False)
        self.early_exit = False
        self.close()    
    
    # need to get rid of the output if exiting early to avoid problems 
    # with the file exiting but not be completed for other tables and 
    # analysis step
    def closeEvent(self, event):
        if self.early_exit:
            os.remove(self.outfile)
        event.accept()
        
   
# we need to do some error checking here.  specifically we need a certain number of
# enrichment values and unique times to actually fit later.
class subject(object):
    def __init__(self, name, time, enrichment):
        self.name = name
        self.times = [time]
        self.enrichments = [enrichment]
    def add_data(self, time, enrichment):
        self.times.append(time)
        self.enrichments.append(enrichment)
        
    #  need to do a small number of basic tests
    def error_checking(self, minimum_time_points):
        if len(set(self.times))< minimum_time_points:
            return "Subject {} has only {} unique time points. {} are requied.".format(
                self.name, len(set(self.times)), minimum_time_points)
        # for now disallow all enrichments for a sample being 0
        if len(set(self.enrichments)) == 1  and self.enrichments[0] == 0.0:
            return "Subject {} has no enrichment.".format(self.name)
        return 1

# test the main
if __name__ == '__main__':
    # needed for windows multiprocessing which will happen at some point
    import sys
    app = QtWidgets.QApplication(sys.argv)
    output_file = ""
    gui_object = EnrichmentWindow(None, ["A", "B", "C"], 3, 5, output_file)
    gui_object.show()
    app.exec_()
