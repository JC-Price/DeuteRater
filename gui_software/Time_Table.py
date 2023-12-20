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
this governs the table the user fills in to provide the time and subject associated with a given mzml file
"""
import os
import csv

from pathlib import Path

from PyQt5 import uic, QtWidgets, QtCore, QtGui

import gui_software.general_table_functions as gtf

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Time_Table.ui")

loaded_ui = uic.loadUiType(ui_file)[0]

#  it is tricky to actually get the header out of the qtablewidget and they
#  need different checks anyway so we'll just declare it here
current_columns = ["Filename", "Time", "Subject ID"]

class TimeWindow(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, filenames = [], outfile = None):
        super(TimeWindow, self).__init__(parent)
        # allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        
        self.setupUi(self)
        
        self.outfile = outfile
        self.filepaths = [Path(f) for f in filenames]
        gtf.set_up_table(self.TimeTable, self.filepaths, rows_are_paths = True)
        
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)
        self.setWindowTitle("Fill in Time Values")
        # make some easy shortcuts. setting up undo and redo are the hardest
        # conly do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(lambda: gtf.Copy_to_Clipboard(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(lambda:gtf.Paste(
                                self.TimeTable))
        
    # we can just write out a .csv or .tsv with the csv module
    # first we need to do some checks
    def check_and_save(self):
        results =[current_columns]
        for r in range(self.TimeTable.rowCount()):
            # filename is not editable so can be ignored
            #current_row = [str(self.TimeTable.item(r,0).text())]
            current_row = [self.filepaths[r].resolve()]
            for i in range(1, len(current_columns)-1):
                test_value = gtf.basic_number_error_check(
                    self.TimeTable.item(r,i).text(),
                    current_columns[i], r)
                if type(test_value) == str:
                    QtWidgets.QMessageBox.information(self, "Error", test_value)
                    return 
                current_row.append(test_value)
            #get the sample group name (no need for it to be )
            test_value, error_code = TimeWindow._basic_string_check(
                    self.TimeTable.item(r,len(current_columns)-1).text(),
                    current_columns[-1], r)
            if error_code:
                QtWidgets.QMessageBox.information(self, "Error", test_value)
                return
            current_row.append(test_value)
            results.append(current_row)
        #  if we've gotten here, we're good to write out
        with open(self.outfile, "w", newline ='') as temp_out:
            writer = csv.writer(temp_out, delimiter ='\t')
            writer.writerows(results)
            #for row in results:
            #    writer.writerow(row)
        self.close()
     
    # check for text errors, like blanks
    @staticmethod    
    def _basic_string_check(text_value, column_name, row_number):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number)
        if text_value == "":
            return "Blank present" + append_to_error,True
        else: 
            return text_value, False
    
# main for testing
if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    output_file = ""
    gui_object = TimeWindow(None, ["A", "B", "C"], output_file)
    gui_object.show()
    app.exec_()
