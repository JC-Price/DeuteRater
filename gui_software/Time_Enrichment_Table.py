# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 08:21:34 2020

@author: JCPrice

The point of this is to hold the time and enrichment table
any adjustments to the table should be made here or in the qt designer
"""
import sys
import os
import csv
import pandas as pd

from PyQt5 import uic, QtWidgets, QtCore, QtGui
import gui_software.general_table_functions as gtf

# when compiling/building for an executable, set all of these to True, otherwise leave as False
# copy "exe_mode = False" and search using ctrl+shift+f to find each instance
exe_mode = False
if exe_mode:
    location = os.path.dirname(os.path.abspath(sys.executable))
else:
    location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Time_Enrichment_Table.ui")

loaded_ui = uic.loadUiType(ui_file)[0]

#  it is tricky to actually get the header out of the qtablewidget and they
#  need different checks anyway, so we'll just declare it here
current_columns = ["Filename", "Time", "Enrichment", "Sample_Group", "Biological_Replicate", "Calculate_N_Value"]


class TimeEnrichmentWindow(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent=None, filenames=[], outfile=None, max_enrichment=None):
        super(TimeEnrichmentWindow, self).__init__(parent)
        # allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint |
                            QtCore.Qt.WindowMaximizeButtonHint)
        self.setupUi(self)

        self.outfile = outfile

        self.max_enrichment = max_enrichment
        self.originalFileNames = filenames
        self.set_up_table(filenames)
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)

        # make some easy shortcuts. setting up undo and redo are the hardest
        # only do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(self.Clear_Contents)
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(self.Clear_Contents)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(self.Copy_to_Clipboard)
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(self.Paste)

        self.setWindowTitle("Time Enrichment Table")

        # resizes columns to fit data
        self.TimeEnrichmentTable.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.TimeEnrichmentTable.resizeColumnsToContents()

    #  initial table set up
    def set_up_table(self, row_names):
        # block signals so if we use undo and redo the basic table is unaffected
        # self.TimeEnrichmentTable.blockSignals(True)
        self.TimeEnrichmentTable.setColumnCount(len(current_columns))
        # self.TimeEnrichmentTable.setHorizontalHeaderLabels(current_columns)
        self.TimeEnrichmentTable.setRowCount(len(row_names))
        for n in range(len(row_names)):
            # only grab the base file name to display on the table. Tried using the following code to make it right justified
            # to see the file name, but it didn't work. - Ben D
            # name_item.setTextAlignment(QtCore.Qt.AlignRight)
            base_name = os.path.basename(row_names[n])
            # make an item
            name_item = QtWidgets.QTableWidgetItem(base_name)
            #  this makes the item uneditable so the user can't mess it up
            #  from https://stackoverflow.com/questions/17104413/pyqt4-how-to-
            #  select-table-rows-and-disable-editing-cells anser 2 accessed 9/25/20
            name_item.setFlags(QtCore.Qt.ItemIsSelectable |
                               QtCore.Qt.ItemIsEnabled)
            self.TimeEnrichmentTable.setItem(n, 0, name_item)
            for i in range(1, self.TimeEnrichmentTable.columnCount()):
                # makes an empty string in all other cells
                self.TimeEnrichmentTable.setItem(n, i, QtWidgets.QTableWidgetItem(""))
        # self.TimeEnrichmentTable.blockSignals(False)

    # we can just write out a .csv or .tsv with the csv module
    # first we need to do some checks
    def check_and_save(self):
        # replace shortened filenames with original full paths - Ben D
        for i in range(len(self.originalFileNames)):
            self.TimeEnrichmentTable.setItem(i, 0, QtWidgets.QTableWidgetItem(self.originalFileNames[i]))

        results = [current_columns]
        for r in range(self.TimeEnrichmentTable.rowCount()):
            # filename is not editable so can be ignored
            current_row = [str(self.TimeEnrichmentTable.item(r, 0).text())]
            for i in range(1, current_columns.index("Sample_Group")):
                test_value = TimeEnrichmentWindow._basic_number_error_check(
                    self.TimeEnrichmentTable.item(r, i).text(),
                    current_columns[i], r)
                if type(test_value) is str:
                    QtWidgets.QMessageBox.information(self, "Error", test_value)
                    return
                if test_value == "inf":
                    from numpy import inf
                    test_value = inf
                # need to block enrichment over 100%
                enrichment_string = self.TimeEnrichmentTable.item(r, 2).text()
                enrichment_value = gtf.basic_number_error_check(enrichment_string, current_columns[2], r, self.max_enrichment, True)
                for value in [test_value, enrichment_value]:
                    if type(value) is str:
                        QtWidgets.QMessageBox.information(self, "Error", value)
                        return
                current_row.append(test_value)
            # get the sample group name and bio rep name (no need for it to be )
            for i in range(current_columns.index("Sample_Group"), len(current_columns)):
                test_value, error_code = TimeEnrichmentWindow._basic_string_check(
                    self.TimeEnrichmentTable.item(r, i).text(),
                    current_columns[i], r)
                if error_code:
                    QtWidgets.QMessageBox.information(self, "Error", test_value)
                    return
                current_row.append(test_value)

            if current_row[i].lower() not in ['yes', 'no']:
                QtWidgets.QMessageBox.information(self, "Reminder",
                                                  "Must answer yes or no for each replicate in the Calculate N-Value column")
                return

            results.append(current_row)
        #  if we've gotten here, we're good to write out
        with open(self.outfile, "w", newline='') as temp_out:
            writer = csv.writer(temp_out, delimiter='\t')
            writer.writerows(results)
            # for row in results:
            #    writer.writerow(row)
        self.close()

    #  everything down to the staticmethods are for table manipulation
    #  they can be shifted to a separate python file if need be
    # if so just replace self.table with a passed in table object

    #  gets the top left selected cell or all selected rows or columns
    def HighlightedCells(self, all_data=False, no_names=True):
        rows = []
        columns = []
        for selected_cell in self.TimeEnrichmentTable.selectedIndexes():
            rows.append(int(selected_cell.row()))
            columns.append(int(selected_cell.column()))
        #  don't let the user paste over or delete the names (copy is fine)
        if no_names:
            columns = [c for c in columns if c != 0]
        if all_data:
            # only get unique rows and columns and sort them (set is unsorted)
            rows = sorted(list(set(rows)))
            columns = sorted(list(set(columns)))
            return rows, columns
        else:
            return rows[0], columns[0]

    def Clear_Contents(self):
        selected_rows, selected_columns = self.HighlightedCells(True)
        for r in selected_rows:
            for c in selected_columns:
                self.TimeEnrichmentTable.item(r, c).setText("")

    def Copy_to_Clipboard(self):
        rows, columns = self.HighlightedCells(all_data=True, no_names=False)
        string_for_clipboard = ""
        for r in rows:
            for c in columns:
                string_for_clipboard += str(self.TimeEnrichmentTable.item(
                    r, c).text())
                string_for_clipboard += "\t"
            #  swap out last tab for a newline
            string_for_clipboard = string_for_clipboard[:-1] + "\n"
        clipboard = QtWidgets.QApplication.clipboard()
        #  trim off the final new line and get it on the clipboard
        clipboard.setText(string_for_clipboard[:-1])

    def Paste(self):
        start_row, start_column = self.HighlightedCells(no_names=False)
        #  do not paste over the first column
        if start_column == 0:
            return
        text = str(QtWidgets.QApplication.clipboard().text())
        if text == "":
            return
        text_lines = text.split("\n")
        text_grid = [line.split("\t") for line in text_lines]
        # excel sometimes adds whitespace to end of clipboard.  must remove
        if text_grid[-1] == ['']:
            del text_grid[-1]
        line_lengths = [len(text_list) for text_list in text_grid]
        needed_columns = max(line_lengths)
        needed_rows = len(text_grid)
        #  check that we don't need more rows or columns that we have
        #  separate lines for easy readability
        if start_row + needed_rows > self.TimeEnrichmentTable.rowCount():
            return
        if start_column + needed_columns > self.TimeEnrichmentTable.columnCount():
            return
        for c in range(needed_columns):
            for r in range(needed_rows):
                self.TimeEnrichmentTable.item(start_row + r, start_column + c).setText(text_grid[r][c])
        self.TimeEnrichmentTable.resizeColumnsToContents()

    #  does some basic error checking for numerical data
    @staticmethod
    def _basic_number_error_check(text_value, column_name, row_number):
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
        # we could also check if it is under some maximum, but for now
        # trust the user
        return num_value

    @staticmethod
    def _basic_string_check(text_value, column_name, row_number):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number)
        if text_value == "":
            return "Blank present" + append_to_error, True
        else:
            return text_value, False


if __name__ == '__main__':
    # needed for windows multiprocessing which will happen at some point
    import sys

    app = QtWidgets.QApplication(sys.argv)
    output_file = "C:\\Software\\Testing\\DeuteRater_Initial\\test_table_out.csv"
    gui_object = TimeEnrichmentWindow(None, ["A", "B", "C"], output_file)
    gui_object.show()
    app.exec_()
