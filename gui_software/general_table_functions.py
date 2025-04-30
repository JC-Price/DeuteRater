# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, Benjamin Driggs,
    Coleman Nielsen, J.C. Price, and Brigham Young University
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

this has various functions for tables in guis, such as making the tables,
finding which cells are highlighted and keyboard shortcuts

"""

from PyQt5 import QtWidgets, QtCore


#  sets up a table, including making some cells uneditable to prevent issues
def set_up_table(table_object, row_names, repeats=1, rows_are_paths=False):
    table_object.setRowCount(len(row_names) * repeats)
    current_row = 0
    for name in row_names:
        for i in range(repeats):
            if rows_are_paths:
                name_item = QtWidgets.QTableWidgetItem(str(name.stem))
            else:
                name_item = QtWidgets.QTableWidgetItem(str(name))
            #  this makes the item uneditable so the user can't mess it up
            # from https://stackoverflow.com/questions/17104413/pyqt4-how-to-
            # select-table-rows-and-disable-editing-cells anser 2 accessed 9/25/20
            name_item.setFlags(QtCore.Qt.ItemIsSelectable |
                               QtCore.Qt.ItemIsEnabled)
            table_object.setItem(current_row, 0, name_item)
            for i in range(1, table_object.columnCount()):
                # makes an empty string in all other cells
                table_object.setItem(current_row, i,
                                     QtWidgets.QTableWidgetItem(""))
            current_row += 1


#  gets the top left selected cell or all selected rows or columns
def HighlightedCells(table_object, all_data=False, no_names=True):
    rows = []
    columns = []
    for selected_cell in table_object.selectedIndexes():
        rows.append(int(selected_cell.row()))
        columns.append(int(selected_cell.column()))
    #  don't let the user paste over or delete the names (copy is fine)
    if no_names:
        columns = [c for c in columns if c != 0]
    if all_data:
        # only get unique rows and colums and sort them (set is unsorted)
        rows = sorted(list(set(rows)))
        columns = sorted(list(set(columns)))
        return rows, columns
    else:
        return rows[0], columns[0]


# clears highlighted cells of their contents
def Clear_Contents(table_object):
    selected_rows, selected_columns = HighlightedCells(table_object, True)
    for r in selected_rows:
        for c in selected_columns:
            table_object.item(r, c).setText("")


# copy highlighted cell values onto the clipboard
def Copy_to_Clipboard(table_object):
    rows, columns = HighlightedCells(table_object,
                                     all_data=True, no_names=False)
    string_for_clipboard = ""
    for r in rows:
        for c in columns:
            string_for_clipboard += str(table_object.item(r, c).text())
            string_for_clipboard += "\t"
        #  swap out last tab for a newline
        string_for_clipboard = string_for_clipboard[:-1] + "\n"
    clipboard = QtWidgets.QApplication.clipboard()
    #  trim off the final new line and get it on the clipboard
    clipboard.setText(string_for_clipboard[:-1])


# paste clipboard into the table. pastes nothing if the thing to be pasted is too large to fit
def Paste(table_object):
    start_row, start_column = HighlightedCells(table_object, no_names=False)
    #  do not paste over the first column
    if start_column == 0:
        return
    text = str(QtWidgets.QApplication.clipboard().text())
    text_lines = text.split("\n")
    text_grid = [line.split("\t") for line in text_lines]
    # excel sometimes adds whitespace to end of clipboard. Must remove
    if text_grid[-1] == ['']: del text_grid[-1]
    line_lengths = [len(text_list) for text_list in text_grid]
    needed_columns = max(line_lengths)
    needed_rows = len(text_grid)
    #  check that we don't need more rows or columns that we have
    #  separate lines for easy readability
    if start_row + needed_rows > table_object.rowCount():
        return
    if start_column + needed_columns > table_object.columnCount():
        return
    for c in range(needed_columns):
        for r in range(needed_rows):
            table_object.item(start_row + r,
                              start_column + c).setText(text_grid[r][c])


# does some basic error checking for numerical data           
def basic_number_error_check(text_value, column_name, row_number, allowed_max=None, check_max=False):
    append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
        column_name, row_number + 1)
    if text_value == "":
        return "Blank present" + append_to_error
    try:
        num_value = float(text_value)
    except ValueError:
        return "Non-numerical value" + append_to_error
    if num_value < 0:
        return "Negative value" + append_to_error

    if check_max:
        if num_value > allowed_max:
            return "Value {} higher than the max allowed error of {}".format(num_value, allowed_max) + append_to_error
    return num_value
