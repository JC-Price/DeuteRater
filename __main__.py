# # -*- coding: utf-8 -*-
# """
# Copyright (c) 2016-2020 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
# All rights reserved.
# Redistribution and use in source and binary forms,
# with or without modification, are permitted provided
# that the following conditions are met:
#     * Redistributions of source code must retain the
#       above copyright notice, this list of conditions
#       and the following disclaimer.
#     * Redistributions in binary form must reproduce
#       the above copyright notice, this list of conditions
#       and the following disclaimer in the documentation
#       and/or other materials provided with the distribution.
#     * Neither the author nor the names of any contributors
#       may be used to endorse or promote products derived
#       from this software without specific prior written
#       permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# """
#
# #$we will of course need to expand things later, but we'll sort that out later
# import os
# import multiprocessing as mp
# import csv
# import pandas as pd
#
# from PyQt5 import uic, QtWidgets, QtCore, QtGui
# from shutil import copyfile, rmtree
#
# from deuteconvert.peaks85 import Peaks85
# from deuteconvert.peaksXplus import PeaksXplus
# from deuteconvert.peaksXpro import PeaksXpro
# from deuterater.extractor import Extractor
# from gui_software.Time_Enrichment_Table import TimeEnrichmentWindow
# from deuterater.theory_preparer import TheoryPreparer
# from deuterater.fraction_new_calculator import FractionNewCalculator
# from deuterater.rate_calculator import RateCalculator
# from utils.useful_classes import deuterater_step, deuteconvert_peaks_required_headers
# import deuterater.settings as settings
# import gui_software.Rate_Settings as rate_settings
#
# #location = os.path.dirname(os.path.abspath(sys.executable))
# location = os.path.dirname(os.path.abspath(__file__))
#
# rate_settings_file = os.path.join(location, "resources","temp_settings.yaml")
# default_rate_settings = os.path.join(location, "resources","settings.yaml")
# id_settings_file = os.path.join(location, "resources",
#                                    "temp_id_settings.yaml")
# default_id_settings = os.path.join(location, "resources",
#                                       "id_settings.yaml")
#
# #$make some basic classes to hold some data.  If need to adjust
# #$output names or columns required from input, do it here
# #$if we add other id types just yank id columns out and make a variable list
# Extract_object = deuterater_step("", ['Precursor Retention Time (sec)',
#                     'Precursor m/z','Identification Charge', 'Sequence',
#                     'Protein ID', "cf"])
# Time_Enrich_object = deuterater_step("time_enrichment_data.tsv",
#     ["Precursor Retention Time (sec)", "Protein ID", "Protein Name", "Precursor m/z",
#     "Identification Charge", "Homologous Proteins", "n_isos",
#     "Sequence", "cf", "abundances", "mzs"
#     ])
# Theory_object = deuterater_step("theory_output.tsv",
#     ["Filename", "Time","Enrichment", "Sample_Group"])
# Fracnew_object = deuterater_step("frac_new_output.tsv", [
#     "Precursor Retention Time (sec)", "Protein ID", "Protein Name", "Precursor m/z",
#     "Identification Charge", "Homologous Proteins", "n_isos", "literature_n",
#     "Sequence", "cf", "abundances", "mzs", "time", "enrichment", "sample_group"])
# #$rate needs reassignment based on settings so we'll read in and make later
# Rate_object = deuterater_step("calculated_rates.csv", [])
#
# step_object_dict = {
#     "Extract":Extract_object,
#     "Provide Time and Enrichment":Time_Enrich_object,
#     "Theory Generation":Theory_object,
#     "Fraction New Calculation":Fracnew_object,
#     "Rate Calculation":Rate_object
#     }
#
# convert_options = {
#     "Peaks 8.5 - Peptides": Peaks85,
#     "Peaks X+ - Peptides": PeaksXplus,
#     "Peaks XPro - Peptides": PeaksXpro,
#     "Template": ""
#     }
#
# convert_needed_headers = {
#     "Peaks 8.5 - Peptides": deuteconvert_peaks_required_headers(
#         ['Accession', '#Peptides', '#Unique', 'Description'],
#         ['Peptide', 'Quality', 'Avg. ppm', 'Start', 'End', 'PTM'],
#         ['Peptide', 'RT range', 'RT mean', 'm/z', 'z', 'Accession', 'PTM']
#         ),
#     "Peaks X+ - Peptides": deuteconvert_peaks_required_headers(
#         ['Accession', '#Peptides', '#Unique', 'Description'],
#         ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
#         ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
#             'RT', 'm/z', 'z', 'Accession']
#         ),
#     "Peaks XPro - Peptides": deuteconvert_peaks_required_headers(
#         ['Accession', '#Peptides', '#Unique', 'Description'],
#         ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
#         ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
#             'RT', 'm/z', 'z', 'Accession']
#         )
#     }
#
#
# default_converter = "Template"
# converter_header = PeaksXplus.correct_header_order
#
# main_file_ui_location = os.path.join(location, "ui_files", "Main_Menu.ui")
# loaded_ui = uic.loadUiType(main_file_ui_location)[0]
# class MainGuiObject(QtWidgets.QMainWindow, loaded_ui):
#     def __init__(self, parent = None):
#         QtWidgets.QMainWindow.__init__(self,parent)
#         #$allows a maximize button
#         self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint |
#                             QtCore.Qt.WindowMaximizeButtonHint)
#
#         self.setupUi(self)
#         make_temp_file(default_rate_settings, rate_settings_file)
#         make_temp_file(default_id_settings, id_settings_file)
#         self.file_loc = location
#
#         #$set up the id file options
#         self.id_file_options.addItems(convert_options.keys())
#         #$set the default value for the converter
#         index = self.id_file_options.findText(default_converter)
#         self.id_file_options.setCurrentIndex(index)
#
#         self.IDFileButton.clicked.connect(self.create_id_file)
#         self.RateCalculationButton.clicked.connect(self.run_rate_workflow)
#         self.actionSettings.triggered.connect(self.change_settings)
#
#         #$make the logo show up
#         #$use of command from http://stackoverflow.com/questions/8687723/pyqthow-do-i-display-a-image-properly
#         #$first answer accesed 5/27/2016
#         myPixmap = QtGui.QPixmap(os.path.join(location, "resources", "Logo.JPG"))
#         self.Logo.setPixmap(myPixmap)
#         self.Logo.setScaledContents(True)
#
#     def Peaks_File_Collection(self, header_checker_object):
#          #$ get the files we need
#         #TODO switch to reading from a folder instead of individual files?$
#         QtWidgets.QMessageBox.information(self, "Info", ("Please select the "
#             "files you would like to turn into a ID file. The order is \""
#             "proteins.csv\", \"protein-peptides.csv\", \"feature.csv\""))
#
#         protein_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
#                 "Choose protein file to Load", self.file_loc, "CSV (*.csv)",
#                 options=QtWidgets.QFileDialog.DontUseNativeDialog)
#         if protein_file == "":
#             return ""
#         else:
#             has_needed_columns = header_checker_object.protein_file_check(protein_file)
#             if not has_needed_columns:
#                 QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
#                     "missing needed columns. Please correct and try again".format(
#                     protein_file)))
#                 return ""
#             self.file_loc = os.path.dirname(protein_file)
#
#         protein_peptide_file, file_type = QtWidgets.QFileDialog.getOpenFileName(
#             self, "Choose protein_peptide file to Load", self.file_loc,
#             "CSV (*.csv)", options=QtWidgets.QFileDialog.DontUseNativeDialog)
#         if protein_peptide_file == "":
#             return ""
#         else:
#             has_needed_columns = header_checker_object.protein_peptide_check(protein_peptide_file)
#             if not has_needed_columns:
#                 QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
#                     "missing needed columns. Please correct and try again".format(
#                     protein_peptide_file)))
#                 return ""
#             self.file_loc = os.path.dirname(protein_peptide_file)
#
#         feature_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
#                 "Choose features file to Load", self.file_loc, "CSV (*.csv)",
#                 options=QtWidgets.QFileDialog.DontUseNativeDialog)
#         if feature_file == "": return ""
#         else:
#             has_needed_columns = header_checker_object.features_check(feature_file)
#             if not has_needed_columns:
#                 QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
#                     "missing needed columns. Please correct and try again".format(
#                     feature_file)))
#                 return ""
#             self.file_loc = os.path.dirname(feature_file)
#         return [protein_file, protein_peptide_file, feature_file]
#
#
#     #$this is to govern the different id file functions
#     def create_id_file(self):
#         id_file_type = str(self.id_file_options.currentText())
#         #$collect any id files needed
#         #$template doesn't need one since it just needs one output
#         if id_file_type in ["Peaks 8.5 - Peptides", "Peaks X+ - Peptides", "Peaks XPro - Peptides"]:
#              input_files = self.Peaks_File_Collection(convert_needed_headers[id_file_type])
#
#         #$id_file_type has to be first or input_files may not be defined
#         if id_file_type != "Template" and input_files =="":
#             return
#         if id_file_type != "Template":
#              #$do the actual calculations
#             converter = convert_options[id_file_type](input_files,
#                 id_settings_file)
#             converter.convert()
#
#         #$get output file
#         QtWidgets.QMessageBox.information(self, "Info", ("Your id file was "
#             "created. Please select the output file location"))
#         while (True):
#             save_file, filetype = QtWidgets.QFileDialog.getSaveFileName(self,
#                 "Provide Save File", self.file_loc, "CSV (*.csv)")
#             if save_file == "": return
#             try:
#                 if id_file_type != "Template":
#                     converter.write(save_file)
#                 else:
#                     df = pd.DataFrame(columns =converter_header )
#                     df.to_csv(save_file, sep=',', index=False)
#                 break
#             except IOError:
#                 QtWidgets.QMessageBox.information(self, "Error",
#                     ("File {} is open in another program. Please close it and "
#                     "try again or select a different file".format(save_file)))
#         self.file_loc = os.path.dirname(save_file)
#         QtWidgets.QMessageBox.information(self, "Success",
#                 "ID file successfully saved")
#
#     def run_rate_workflow(self):
#         #$will need some settings
#         settings.load(rate_settings_file)
#
#         #$first we need to check which steps are checked
#         worklist = self.check_table_checklist()
#         #$ only proceed if we have a
#         if worklist == []:
#             QtWidgets.QMessageBox.information(self, "Error", ("No options were "
#                 "checked. Please check steps to performed and try again"))
#             return
#         elif type(worklist) == str:
#             QtWidgets.QMessageBox.information(self, "Error", worklist)
#             return
#         #$second we need to an output folder and check it for output folder
#         QtWidgets.QMessageBox.information(self, "Info", ("Please select folder "
#                 "for output"))
#         output_folder = QtWidgets.QFileDialog.getExistingDirectory(
#             self,
#             "Select an Output Folder",
#             self.file_loc,
#             QtWidgets.QFileDialog.ShowDirsOnly)
#         if output_folder == "": return
#         #$change location we start asking for things at
#         #$don't change since all output is going in here
#         self.file_loc = output_folder
#         self.make_folder(output_folder, non_graph=True)
#
#         #$then need to check if the files exist. if so warn the user. function
#         no_extract_list = [w for w in worklist if w !="Extract"]
#         outputs_to_check = []
#         for worklist_step in no_extract_list:
#             step_object_dict[worklist_step].complete_filename(self.file_loc)
#             outputs_to_check.append(step_object_dict[worklist_step].full_filename)
#         #$if should only fail if an extract only, but that may occur
#         if outputs_to_check != []:
#             proceed = self.check_file_removal(outputs_to_check)
#             if not proceed:
#                 return
#
#         #$now we need to get input and do the work. each step can only occur
#         #$once and they occur in order. so we will write them in order
#         #todo$ see if we can compress the code and make sure it is readable
#         previous_output_file = ""
#         extracted_files =[]
#         make_table_in_order = True
#         for analysis_step in worklist:
#             if analysis_step == "Extract":
#                 #$no if for this one, if extract is here it is the start
#                 id_file = self.collect_single_file("ID", "Extract", "CSV (*.csv)")
#                 if id_file == "": return
#                 #$always check if is good since it is first
#                 infile_is_good = self.check_input(step_object_dict["Extract"],
#                                         id_file)
#                 if not infile_is_good:  return
#
#                 mzml_files = self.collect_multiple_files("Centroided Data",
#                                 "Extract", "mzML (*.mzML)")
#                 if mzml_files ==[]: return
#
#                 mzml_filenames = [os.path.basename(filename) for filename in
#                                  mzml_files]
#                 extracted_files = [filename.replace(".mzML", ".tsv") for
#                                         filename in mzml_filenames]
#                 extracted_files = [os.path.join(output_folder, filename) for
#                                    filename in extracted_files]
#
#                 proceed = self.check_file_removal(extracted_files)
#                 if not proceed:  return
#                 #$need to run the table if necessary. taken from the
#                 #$"Provide Time and Enrichment" elif
#                 if "Provide Time and Enrichment" in worklist:
#                     previous_output_file = step_object_dict[
#                         "Provide Time and Enrichment"].full_filename
#                     self.get_data_table = TimeEnrichmentWindow(self,
#                             extracted_files, previous_output_file)
#                     self.get_data_table.exec_()
#                     #$don't make the table twice
#                     make_table_in_order = False
#                     #$now that the table is done we need to confirm the user
#                     #$hit the proceed button on the table (same check as in
#                     #$elif analysis_step == "Theory Generation" )
#                     if not os.path.exists(previous_output_file): return
#                 #$ modified from the extract-dir argument from the command line
#                 for m in range(len(mzml_files)):
#                     extractor = Extractor(
#                         id_path = os.path.join(self.file_loc, id_file),
#                         mzml_path = mzml_files[m],
#                         out_path = extracted_files[m],
#                         settings_path = rate_settings_file
#                         )
#                     extractor.load()
#                     extractor.run()
#                     extractor.write()
#             elif analysis_step == "Provide Time and Enrichment" and make_table_in_order:
#                 #$if coming right after a list
#                 if extracted_files == []:
#                     extracted_files = self.collect_multiple_files(
#                         "Extracted Data",
#                         "Provide Time and Enrichment",
#                         "TSV (*.tsv)"
#                         )
#                     if extracted_files == []: return
#                     #$ensure the input files are good. only need to deal with
#                     #$if the user just selected
#                     for e_file in extracted_files:
#                         infile_is_good = self.check_input(
#                             step_object_dict["Provide Time and Enrichment"],
#                             e_file)
#                         if not infile_is_good: return
#
#                 #$ now that we have the extracted files we can make a table
#                 #$the talbe will handle the output
#                 previous_output_file = step_object_dict[
#                     "Provide Time and Enrichment"].full_filename
#                 self.get_data_table = TimeEnrichmentWindow(self,
#                         extracted_files, previous_output_file)
#                 self.get_data_table.exec_()
#             elif analysis_step == "Theory Generation":
#                 #$since the files are in the table can just read that in
#                 if previous_output_file == "":
#                     previous_output_file = self.collect_single_file(
#                         "time and enrichment",
#                         "Theory Generation",
#                         "spreadsheet (*.csv *.tsv)"
#                     )
#                     if previous_output_file == "": return
#                     infile_is_good = self.check_input(
#                         step_object_dict["Theory Generation"],
#                         previous_output_file)
#                     if not infile_is_good: return
#                 #$else is to deal with a failed write from the previous table
#                 #$ don't need an error message just return
#                 elif not os.path.exists(previous_output_file): return
#
#                 #$final check to see if all of the files in the input table
#                 #$still exist.  don't want to error out in the middle of
#                 #$multiprocessing
#                 final_proceed =  self.check_files_from_files(
#                     previous_output_file, 0)
#                 if not final_proceed: return
#
#                 theorist = TheoryPreparer(
#                     enrichment_path=previous_output_file,
#                     out_path= step_object_dict["Theory Generation"].full_filename,
#                     settings_path=rate_settings_file
#                 )
#                 theorist.prepare()
#                 theorist.write()
#                 previous_output_file = step_object_dict["Theory Generation"].full_filename
#             elif analysis_step == "Fraction New Calculation":
#                 if previous_output_file == "":
#                     previous_output_file = self.collect_single_file(
#                         "theoretical output",
#                         "Fraction New Calculation",
#                         "spreadsheet (*.csv *.tsv)"
#                     )
#                     if previous_output_file == "": return
#
#                     infile_is_good = self.check_input(
#                         step_object_dict["Fraction New Calculation"],
#                         previous_output_file)
#                     if not infile_is_good: return
#                 #$not sure why this would happen but we'll put it here
#                 #$to avoid future error
#                 elif not os.path.exists(previous_output_file): return
#                 fnewcalc = FractionNewCalculator(
#                     model_path=previous_output_file,
#                     out_path=step_object_dict["Fraction New Calculation"].full_filename,
#                     settings_path=rate_settings_file
#                 )
#                 fnewcalc.generate()
#                 if fnewcalc.error != "":
#                     QtWidgets.QMessageBox.information(self, "Error", fnewcalc.error)
#                     return
#                 fnewcalc.write()
#                 previous_output_file = step_object_dict[
#                     "Fraction New Calculation"].full_filename
#             elif analysis_step == "Rate Calculation":
#                 if previous_output_file == "":
#                     previous_output_file = self.collect_single_file(
#                         "fraction new",
#                         "Rate Calculation",
#                         "spreadsheet (*.csv *.tsv)"
#                     )
#                     if previous_output_file == "": return
#                     #$need to ensure that we have proper colmns which varies
#                     #$by setting
#                     needed_columns = [
#                         settings.peptide_analyte_id_column,
#                         settings.peptide_analyte_name_column,
#                         "sample_group"]
#                     if settings.use_abundance:
#                         needed_columns.extend(["afn", "frac_new_abunds_std_dev"])
#                     if settings.use_neutromer_spacing:
#                         needed_columns.extend(["nsfn", "frac_new_mzs_std_dev"])
#                     if settings.use_abundance and settings.use_neutromer_spacing:
#                         needed_columns.extend(["cfn", "frac_new_combined_std_dev"])
#                     step_object_dict["Rate Calculation"].required_columns =needed_columns
#                     infile_is_good = self.check_input(
#                         step_object_dict["Rate Calculation"],
#                         previous_output_file)
#                     if not infile_is_good: return
#
#                 #$need to get a graph folder and ensure it exists
#                 #$don't worry about overwriting files
#                 GraphFolder = os.path.join(self.file_loc, "Graph_Folder")
#                 self.make_folder(GraphFolder)
#                 ratecalc = RateCalculator(
#                     model_path = previous_output_file,
#                     out_path = step_object_dict["Rate Calculation"].full_filename,
#                     graph_folder = GraphFolder,
#                     settings_path = rate_settings_file
#                 )
#                 ratecalc.calculate()
#                 ratecalc.write()
#         QtWidgets.QMessageBox.information(self, "Success",
#                 "Analysis completed successfully")
#
#
#     def check_table_checklist(self):
#         current_worklist =[]
#         error_check = []
#         for i in range(self.RequestedStepsTable.rowCount()):
#             if self.RequestedStepsTable.item(i,0).checkState() == QtCore.Qt.Checked:
#                 current_worklist.append(str(
#                     self.RequestedStepsTable.item(i,0).text()))
#                 error_check.append(i)
#         #$the point of the error check is to ensure ther are no gaps in the
#         #$checklist since then we will be missing critical info for the next
#         #$step
#         proper_length = max(error_check) - min(error_check) + 1
#         if len (current_worklist) < proper_length:
#             return ("There are gaps in your worklist.  Please check all boxes "
#                     "between the first and last checked box")
#         return current_worklist
#
#     #$quick way to get a response to a question to the user.  exiting is a "No"
#     def question_for_user(self, message):
#         response = QtWidgets.QMessageBox.question(self, "Question", message,
#                         QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
#         if response == QtWidgets.QMessageBox.Yes:
#             return True
#         else:
#             return False
#
#     def collect_single_file(self, to_load, step_name, load_type):
#         QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
#                 "file to load for {} step".format(to_load, step_name)))
#         filename, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
#                 "Choose {} file to load".format(to_load), self.file_loc,
#                 load_type, options=QtWidgets.QFileDialog.DontUseNativeDialog)
#         return filename
#
#     def collect_multiple_files(self, to_load, step_name, load_type):
#         QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
#                 "files to load for {} step".format(to_load, step_name)))
#         filenames, file_type = QtWidgets.QFileDialog.getOpenFileNames(self,
#                 "Choose {} file to load".format(to_load), self.file_loc,
#                 load_type, options=QtWidgets.QFileDialog.DontUseNativeDialog)
#         return filenames
#
#     def change_settings(self):
#         self.set_menu = rate_settings.Rate_Setting_Menu(self, rate_settings_file)
#         self.set_menu.exec()
#
#     #$ we have some cases where we need to remove files that we will create
#     #$later (and so would be overwritten anyway). we'll just do error messages
#     #$ and so on here.
#     def check_file_removal(self, list_of_filenames):
#         files_to_remove = []
#         open_files = []
#         #$find files that exist
#         for filename in list_of_filenames:
#             if os.path.exists(filename):
#                 files_to_remove.append(filename)
#         #$let the user decide if we should continue.
#         if files_to_remove != []:
#             proceed = self.question_for_user(("The following files will be "
#                 "overwritten:\n{}\nDo you still wish to proceed?".format(
#                  ",\n".join(files_to_remove))))
#             if not proceed:
#                  return False
#             for filename in files_to_remove:
#                  try:
#                      os.remove(filename)
#                  except PermissionError:
#                      open_files.append(filename)
#             if open_files != []:
#                 QtWidgets.QMessageBox.information(self, "Error",
#                     ("The following files cannot be overwritten:\n{}\n"
#                      "They are likely open in another program. please close "
#                      "and try again.".format(",\n".join(open_files)))
#                     )
#                 return False
#         #$will return true if no files already exist or the user wants to
#         #$overwrite and they can be removed so we have permission
#         return True
#
#     def check_input (self, relevant_object, filename):
#         has_needed_columns = relevant_object.check_input_file(filename)
#         if not has_needed_columns:
#             QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
#                 "missing needed columns. Please correct and try again".format(
#                 filename)))
#         return has_needed_columns
#
#     #$there are cases (specifically the table going to theory) where there
#     #$is a chance that the files referenced in a id file may not exist
#     #$this will check for that
#     def check_files_from_files(self, input_file, filename_column):
#         with open (input_file, 'r') as infile:
#             if input_file[-4:] == ".tsv":
#                 reader = csv.reader(infile, delimiter = "\t")
#             elif input_file[-4:] == ".csv":
#                 reader = csv.reader(infile)
#             next(reader)#$skip header
#             #$may wish to add a check that rows are of appropriate length
#             for row in reader:
#                 if not os.path.exists(row[filename_column]):
#                     QtWidgets.QMessageBox.information(self, "Error", ("File {} "
#                         "could not be found. Please correct input file and try "
#                         "again".format(row[filename_column]
#                         )))
#                     return False
#         return True
#
#     #$ensures a folder exists by making it if it does not.
#     def make_folder(self, folder, non_graph = False):
#         if non_graph: #$ don't care aobut overwriting graph folder, that is necessary.   the main output folder might have necessary things in it the user wants to keep
#             if not os.path.isdir(folder):
#                 os.makedirs(folder)
#             """
#             #$this is for is we want to delete the entire folder.  this particular use is far too crude.  Any subfolders or files
#             #$even unrelated to DeuteRater will be deleted.  If we need to put this back in we need to target only DeuteRater files
#             #$that should be unnecessary because we have the the overwrite warning
#             if os.path.isdir(folder):
#                 answer = self.question_for_user(f"Would you like to delete the contents of folder {folder}?")
#                 if answer:
#                     rmtree(folder)
#                     os.makedirs(folder)
#             else:
#                 os.makedirs(folder)"""
#         else:
#             if os.path.isdir(folder):
#                 rmtree(folder)
#             os.makedirs(folder)
#
#     @staticmethod
#     def _get_file_names(folder_loc, operations, object_dict):
#         needed_files = []
#         for o in operations:
#             if o != "Extract":
#                 full_filename = os.path.join(folder_loc,
#                                              object_dict.output_filename)
#                 needed_files.append(full_filename)
#         return needed_files
#
# #$since we have to load settings in each file, and need a way to adjust
# #$settings, we'll
# def make_temp_file(filename, new_filename):
#     copyfile(filename, new_filename)
#
# if __name__ == '__main__':
#     #$needed for windows multiprocessing which will happen at some point
#     import sys
#     mp.freeze_support()
#     app = QtWidgets.QApplication(sys.argv)
#     gui_object = MainGuiObject(None)
#     gui_object.show()
#     app.exec_()

if __name__ == '__main__':
    with open("C:\\Users\\benny\\OneDrive\\Desktop\\logs.txt", 'w') as outf:
        outf.write("__main__.py - Running main method in main_gui.py from __main__.py\n")
        from main_gui import main
        main()

