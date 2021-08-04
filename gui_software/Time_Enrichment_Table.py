## -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor, J.C. Price, and Brigham Young University
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
import csv
import pandas as pd

from pathlib import Path


from PyQt5 import uic, QtWidgets, QtCore, QtGui


#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Time_Enrichment_Table.ui")


loaded_ui = uic.loadUiType(ui_file)[0]

#$ it is tricky to actually get the header out of the qtablewidget and they
#$ need different checks anyway so we'll just declare it here
current_columns = ["Filename", "Time", "Enrichment", "Sample_Group"]

enrich_col_loc = current_columns.index("Enrichment")

class TimeEnrichmentWindow(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, filenames = [], outfile = None):
        super(TimeEnrichmentWindow, self).__init__(parent)
        #$allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        self.setupUi(self)
        
        self.outfile = outfile
        self.filepaths = [Path(f) for f in filenames]
        
        self.set_up_table()
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)
        
        #$make some easy shortcuts. setting up undo and redo are the hardest
        #$conly do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(self.Clear_Contents)
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(self.Clear_Contents)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(self.Copy_to_Clipboard)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(self.Paste)
        
    #$ initial table set up
    def set_up_table(self):
        #$block signals so if we use undo and redo the basic table is unaffected
        #self.TimeEnrichmentTable.blockSignals(True)
        self.TimeEnrichmentTable.setRowCount(len(self.filepaths))
        for n in range(len(self.filepaths)):
            #$make an item
            name_item = QtWidgets.QTableWidgetItem(self.filepaths[n].stem)
            #$ this makes the item uneditable so the user can't mess it up
            #$ from https://stackoverflow.com/questions/17104413/pyqt4-how-to-
            #$ select-table-rows-and-disable-editing-cells anser 2 accessed 9/25/20
            name_item.setFlags(QtCore.Qt.ItemIsSelectable | 
                               QtCore.Qt.ItemIsEnabled)
            self.TimeEnrichmentTable.setItem(n, 0, name_item)
            for i in range(1, self.TimeEnrichmentTable.columnCount()):
                #$makes an empty string in all other cells
                self.TimeEnrichmentTable.setItem(n, i, QtWidgets.QTableWidgetItem(""))
        #self.TimeEnrichmentTable.blockSignals(False) 
        
    #$we can just write out a .csv or .tsv with the csv module
    #$first we need to do some checks
    def check_and_save(self):
        results =[current_columns]
        for r in range(self.TimeEnrichmentTable.rowCount()):
            #$filename is not editable so can be ignored
            #current_row = [str(self.TimeEnrichmentTable.item(r,0).text())]
            current_row = [self.filepaths[r].resolve()]
            for i in range(1, len(current_columns)-1):
                if i != enrich_col_loc: 
                    provide_as_percent = False
                else:
                    provide_as_percent = True
                #$r + 1 because it is for an error message.  it will help the user
                test_value = TimeEnrichmentWindow._basic_number_error_check(
                    self.TimeEnrichmentTable.item(r,i).text(),
                    current_columns[i], r+1, provide_as_percent)
                if type(test_value) == str:
                    QtWidgets.QMessageBox.information(self, "Error", test_value)
                    return 
                current_row.append(test_value)
            #get the sample group name (no need for it to be )
            test_value, error_code = TimeEnrichmentWindow._basic_string_check(
                    self.TimeEnrichmentTable.item(r,len(current_columns)-1).text(),
                    current_columns[-1], r+1)
            if error_code:
                QtWidgets.QMessageBox.information(self, "Error", test_value)
                return
            current_row.append(test_value)
            results.append(current_row)
        #$ if we've gotten here, we're good to write out
        with open(self.outfile, "w", newline ='') as temp_out:
            writer = csv.writer(temp_out, delimiter ='\t')
            writer.writerows(results)
            #for row in results:
            #    writer.writerow(row)
        self.close()
     
    #$ everything down to the staticmethods are for table manipulation
    #$ they can be shifted to a separate python file if need be
    #$if so just replace self.table with a passed in table object
    
    
    #$ gets the top left selected cell or all selected rows or columns
    def HighlightedCells(self, all_data = False, no_names = True):
        rows = []
        columns = []
        for selected_cell in self.TimeEnrichmentTable.selectedIndexes():
            rows.append(int(selected_cell.row()))
            columns.append(int(selected_cell.column()))
        #$ don't let the user paste over or delete the names (copy is fine)
        if no_names:
            columns = [c for c in columns if c!= 0]
        if all_data:
            #$only get unique rows and colums and sort them (set is unsorted)
            rows = sorted(list(set(rows)))
            columns = sorted(list(set(columns)))
            return rows, columns
        else:
            return rows[0], columns[0]
        
    def Clear_Contents(self):
        selected_rows, selected_columns = self.HighlightedCells(True)
        for r in selected_rows:
            for c in selected_columns:
                self.TimeEnrichmentTable.item(r,c).setText("")
                
    def Copy_to_Clipboard(self):
        rows, columns = self.HighlightedCells(all_data = True, no_names = False)
        string_for_clipboard = ""
        for r in rows:
            for c in columns:
                string_for_clipboard += str(self.TimeEnrichmentTable.item(
                    r,c).text())
                string_for_clipboard += "\t"
            #$ swap out last tab for a newline
            string_for_clipboard = string_for_clipboard[:-1] + "\n"
        clipboard = QtWidgets.QApplication.clipboard()
        #$ trim off the final new line and get it on the clipboard
        clipboard.setText(string_for_clipboard[:-1]) 
        
    def Paste(self):
        start_row, start_column = self.HighlightedCells(no_names = False)
        #$ do not paste over the first column
        if start_column == 0:
            return
        text = str(QtWidgets.QApplication.clipboard().text())
        text_lines = text.split("\n")
        text_grid = [line.split("\t") for line in text_lines]
        #excel sometimes adds whitespace to end of clipboard.  must remove 
        if text_grid[-1] == ['']: del text_grid[-1]
        line_lengths = [len(text_list) for text_list in text_grid]
        needed_columns = max(line_lengths)
        needed_rows =len(text_grid)
        #$ check that we don't need more rows or columns that we have
        #$ separate lines for easy readability
        if start_row + needed_rows > self.TimeEnrichmentTable.rowCount():
            return
        if start_column + needed_columns > self.TimeEnrichmentTable.columnCount():
            return
        for c in range(needed_columns):
            for r in range(needed_rows):
                self.TimeEnrichmentTable.item(start_row +r,
                    start_column + c).setText(text_grid[r][c])
        
            
        
        
    #$does some basic error checking for numerical data 
    @staticmethod
    def _basic_number_error_check(text_value, column_name, row_number, percent_as_decimal):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number)
        if text_value == "":
            return "Blank present" + append_to_error
        try:
            num_value = float(text_value)
        except ValueError:
            return "Non-numerical value" + append_to_error
        if num_value < 0:
            return "Negative value" + append_to_error
        if percent_as_decimal and num_value >1: #$allow up to 100%
            return "Enrichment is a decimal not a percent. Value is too large" + append_to_error
        
        #$we could also check if it is under some maximum, but for now
        #$trust the user
        return num_value
    @staticmethod    
    def _basic_string_check(text_value, column_name, row_number):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number)
        if text_value == "":
            return "Blank present" + append_to_error,True
        else: 
            return text_value, False
        

if __name__ == '__main__':
    #$needed for windows multiprocessing which will happen at some point
    import sys
    app = QtWidgets.QApplication(sys.argv)
    output_file = "C:\\Software\\Testing\\DeuteRater_Initial\\test_table_out.tsv"
    gui_object = TimeEnrichmentWindow(None, ["C:\\test_folder_name\\A.tsv", "C:\\test_folder_name\\B.tsv", "C:\\test_folder_name\\C.tsv"], output_file)
    gui_object.show()
    app.exec_()
