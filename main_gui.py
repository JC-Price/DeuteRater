# -*- coding: utf-8 -*-
"""
original DeuteRater: Copyright (c) 2016 Bradley Naylor, Michael Porter, J.C. Price, and Brigham Young University

current version created 2021 updates by Bradley Naylor, Kyle Cutler, Chad Quilling, J.C. Price and Brigham Young University
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

# $we will of course need to expand things later, but we'll sort that out later
import os
import multiprocessing as mp
import csv
import pandas as pd
from tqdm import tqdm

from PyQt5 import uic, QtWidgets, QtCore, QtGui
from shutil import copyfile, rmtree

from deuteconvert.peaks85 import Peaks85
from deuteconvert.peaksXplus import PeaksXplus
from deuteconvert.peaksXpro import PeaksXpro
from deuterater.extractor import Extractor
import deuteconvert.peptide_utils as peputils
from gui_software.Time_Enrichment_Table import TimeEnrichmentWindow
from deuterater.theory_preparer import TheoryPreparer
from deuterater.fraction_new_calculator import FractionNewCalculator
from deuterater.rate_calculator import RateCalculator
from utils.useful_classes import deuterater_step, deuteconvert_peaks_required_headers
import deuterater.settings as settings
import gui_software.Rate_Settings as rate_settings
import gui_software.Converter_Settings as guide_settings
import deuteconvert.settings as converter_settings

# location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.abspath(__file__))

rate_settings_file = os.path.join(location, "resources", "temp_settings.yaml")
default_rate_settings = os.path.join(location, "resources", "settings.yaml")
id_settings_file = os.path.join(location, "resources",
                                   "temp_id_settings.yaml")
default_id_settings = os.path.join(location, "resources",
                                      "id_settings.yaml")

# $make some basic classes to hold some data.  If need to adjust
# $output names or columns required from input, do it here
# $if we add other id types just yank id columns out and make a variable list
#$if you change any of these from csv to tsv or vice versa remember to adjust the write function for the relevant class
Extract_object = deuterater_step("", ['Precursor Retention Time (sec)',
                                      'Precursor m/z', 'Identification Charge', 'Sequence',
                                      'Protein ID', "cf"])
Time_Enrich_object = deuterater_step("time_enrichment_data.tsv",
                                     [
                                         "Precursor Retention Time (sec)", "Protein ID", "Protein Name",
                                         "Precursor m/z",
                                         "Identification Charge", "Homologous Proteins", "n_isos", "literature_n",
                                         "Sequence", "cf", "abundances", "mzs"
                                     ])
Theory_object = deuterater_step("theory_output.tsv",
                                ["Filename", "Time", "Enrichment", "Sample_Group"])
Fracnew_object = deuterater_step("frac_new_output.tsv", [
    "Precursor Retention Time (sec)", "Protein ID", "Protein Name", "Precursor m/z",
    "Identification Charge", "Homologous Proteins", "n_isos", "literature_n",
    "Sequence", "cf", "abundances", "mzs", "time", "enrichment", "sample_group"])
# $rate needs reassignment based on settings so we'll read in and make later
Rate_object = deuterater_step("calculated_rates.csv", [])

#$we have an extra file from the rate output. just carve out an exception. if another step gets another output,
#$then adjust the deuterater step object
extra_rate_file = "calculated_rates_datapoints.tsv"

step_object_dict = {
    "Extract": Extract_object,
    "Provide Time and Enrichment": Time_Enrich_object,
    "Theory Generation": Theory_object,
    "Fraction New Calculation": Fracnew_object,
    "Rate Calculation": Rate_object
}

convert_options = {
    "Peaks 8.5 - Peptides": Peaks85,
    "Peaks X+ - Peptides": PeaksXplus,
    "Peaks XPro - Peptides": PeaksXpro,
    "Template": ""
}
convert_needed_headers = {
    "Peaks 8.5 - Peptides": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'Quality', 'Avg. ppm', 'Start', 'End', 'PTM'],
        ['Peptide', 'RT range', 'RT mean', 'm/z', 'z', 'Accession', 'PTM']
        ),
    "Peaks X+ - Peptides": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
        ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
            'RT', 'm/z', 'z', 'Accession']
        ),
    "Peaks XPro - Peptides": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
        ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
            'RT', 'm/z', 'z', 'Accession']
        )
    }

#$columns for the extractor that absolutely need data (others can autofill or are just for user information)
#$most of the columns actually aren't necessary to get a result but will either prevent any extraction from happening (no data output file)
#$or will likely cause an error later in the process.
required_data_extractor_data = ["Sequence", "Protein ID", "Precursor Retention Time (sec)", "Precursor m/z", "Identification Charge"]
autofill_columns = ["Peptide Theoretical Mass", "cf", "literature_n"]
#$need to autofill: Peptide Theoretical Mass, cf? 	neutromers_to_extract, literature_n



default_converter = "Peaks XPro - Peptides"
# TODO$ may need to adjust the header or shove in the n-value calculator
converter_header = PeaksXplus.correct_header_order

main_file_ui_location = os.path.join(location, "ui_files", "Main_Menu.ui")
loaded_ui = uic.loadUiType(main_file_ui_location)[0]


class MainGuiObject(QtWidgets.QMainWindow, loaded_ui):
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        
        # $allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint |
                            QtCore.Qt.WindowMaximizeButtonHint)
        
        self.setupUi(self)
        
        make_temp_file(default_rate_settings, rate_settings_file)
        make_temp_file(default_id_settings, id_settings_file)
        self.file_loc = location
        
        #$set up the id file options
        self.id_file_options.addItems(convert_options.keys())
        # $set the default value for the converter
        index = self.id_file_options.findText(default_converter)
        self.id_file_options.setCurrentIndex(index)
        
        self.IDFileButton.clicked.connect(self.create_id_file)
        self.RateCalculationButton.clicked.connect(self._calc_rates)
        self.actionSettings.triggered.connect(self.change_settings)
        self.actionID_File_Settings.triggered.connect(self.change_converter_settings)
        
        # $make the logo show up
        # $use of command from http://stackoverflow.com/questions/8687723/pyqthow-do-i-display-a-image-properly
        # $first answer accesed 5/27/2016
        myPixmap = QtGui.QPixmap(os.path.join(location, "resources", "Logo.JPG"))
        self.Logo.setPixmap(myPixmap)
        self.Logo.setScaledContents(True)
    
    def Peaks_File_Collection(self, header_checker_object):
        # $ get the files we need
        # TODO switch to reading from a folder instead of individual files?$
        QtWidgets.QMessageBox.information(self, "Info", ("Please select the "
            "files you would like to turn into a ID file. The order is \""
            "proteins.csv\", \"protein-peptides.csv\", \"feature.csv\""))
        
        protein_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
                                                                        "Choose protein file to Load", self.file_loc,
                                                                        "CSV (*.csv)",
                                                                        options=QtWidgets.QFileDialog.DontUseNativeDialog)
        if protein_file == "":
            return ""
        else:
            has_needed_columns = header_checker_object.protein_file_check(protein_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    protein_file)))
                return ""
            self.file_loc = os.path.dirname(protein_file)
        protein_peptide_file, file_type = QtWidgets.QFileDialog.getOpenFileName(
            self, "Choose protein_peptide file to Load", self.file_loc,
            "CSV (*.csv)", options=QtWidgets.QFileDialog.DontUseNativeDialog)
        if protein_peptide_file == "":
            return ""
        else:
            has_needed_columns = header_checker_object.protein_peptide_check(protein_peptide_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    protein_peptide_file)))
                return ""
            self.file_loc = os.path.dirname(protein_peptide_file)
        feature_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
                                                                        "Choose features file to Load", self.file_loc,
                                                                        "CSV (*.csv)",
                                                                        options=QtWidgets.QFileDialog.DontUseNativeDialog)
        if feature_file == "":
            return ""
        else:
            has_needed_columns = header_checker_object.features_check(feature_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    feature_file)))
                return ""
            self.file_loc = os.path.dirname(feature_file)
        return [protein_file, protein_peptide_file, feature_file]
    
    #$this is to govern the different id file functions
    # TODO$ when we have added lipid data adjust template and the function calls
    def create_id_file(self):
        id_file_type = str(self.id_file_options.currentText())
        #$collect any id files needed
        # $template doesn't need one since it just needs one output
        if id_file_type in ["Peaks 8.5 - Peptides", "Peaks X+ - Peptides", "Peaks XPro - Peptides"]:
             input_files = self.Peaks_File_Collection(convert_needed_headers[id_file_type])
        
        #$id_file_type has to be first or input_files may not be defined
        if id_file_type != "Template" and input_files =="":
            return
        if id_file_type != "Template":
            # $do the actual calculations
            converter = convert_options[id_file_type](input_files,
                id_settings_file)
            converter.convert()
        
        # $get output file
        QtWidgets.QMessageBox.information(self, "Info", ("Your id file was "
                                                         "created. Please select the output file location"))
        while (True):
            save_file, filetype = QtWidgets.QFileDialog.getSaveFileName(self,
                                                                        "Provide Save File", self.file_loc,
                                                                        "CSV (*.csv)")
            if save_file == "": return
            try:
                if id_file_type != "Template":
                    converter.write(save_file)
                else:
                    df = pd.DataFrame(columns=converter_header)
                    df.to_csv(save_file, sep=',', index=False)
                break
            except IOError:
                QtWidgets.QMessageBox.information(self, "Error",
                                                  ("File {} is open in another program. Please close it and "
                                                   "try again or select a different file".format(save_file)))
        self.file_loc = os.path.dirname(save_file)
        QtWidgets.QMessageBox.information(self, "Success",
                "ID file successfully saved")
        
    #$adjusts the text fo the "Rate Calculation" button
    def _calc_rates(self):
        try:
            
            self.RateCalculationButton.setText("Currently Processing... Please Wait.")
            self.run_rate_workflow()
        finally:
            self.RateCalculationButton.setText("Rate Calculation")
    
    def run_rate_workflow(self):
        # $will need some settings
        settings.load(rate_settings_file)
        
        # $first we need to check which steps are checked
        worklist = self.check_table_checklist()
        # $ only proceed if we have a
        if worklist == []:
            QtWidgets.QMessageBox.information(self, "Error", ("No options were "
                                                              "checked. Please check steps to performed and try again"))
            return
        elif type(worklist) == str:
            QtWidgets.QMessageBox.information(self, "Error", worklist)
            return
        # $second we need to an output folder and check it for output folder
        QtWidgets.QMessageBox.information(self, "Info", ("Please select folder "
                                                         "for output"))
        output_folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select an Output Folder",
            self.file_loc,
            QtWidgets.QFileDialog.ShowDirsOnly)
        if output_folder == "": return
        # $change location we start asking for things at
        # $don't change since all output is going in here
        self.file_loc = output_folder
        # MainGuiObject._make_folder(output_folder)
        
        #don't care if overwrite rate_settings.yaml but should check if want to use the settings already in the folder
        if os.path.exists(os.path.join(output_folder, "rate_settings.yaml")):
            comp_result = settings.compare(rate_settings_file, os.path.join(output_folder, "rate_settings.yaml"))
            if comp_result != "MATCH":
                if comp_result == "Mismatched Keys":
                    qBox = QtWidgets.QMessageBox(self)
                    qBox.setWindowTitle("Question")
                    question = "A settings file already exists in this output folder. Would you like to use those settings,or overwrite them?"
                    qBox.setText(question)
                    qBox.setIcon(QtWidgets.QMessageBox.Question)
                    qBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)

                    yButton = qBox.button(QtWidgets.QMessageBox.Yes)
                    yButton.setText("Use Settings")
                    nButton = qBox.button(QtWidgets.QMessageBox.No)
                    nButton.setText("Overwrite")
                    response = qBox.exec_()
                    if response == QtWidgets.QMessageBox.Yes:
                        settings.load(os.path.join(output_folder, "rate_settings.yaml"))
                        settings.freeze(rate_settings_file)
                    elif response == QtWidgets.QMessageBox.No:
                        if self.check_file_removal([os.path.join(output_folder, "rate_settings.yaml")], ask_permission = False):
                            settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
                        else:
                            return
                    else:
                        return
                else:
                    #$no point asking if we can delete the file if it is the same anyway.  still want to overwrite as part of checking permissions.
                    #$may not overwrite later
                    if self.check_file_removal([os.path.join(output_folder, "rate_settings.yaml")], ask_permission = False):
                        settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
                    else:
                        return
        else:
            settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
        
        # $then need to check if the files exist. if so warn the user. function
        no_extract_list = [w for w in worklist if w != "Extract"]
        outputs_to_check = []
        #$deal with the extra detail file added for rate calculation
        if "Rate Calculation" in no_extract_list:
            outputs_to_check.append(os.path.join(output_folder, extra_rate_file))
        
        for worklist_step in no_extract_list:
            step_object_dict[worklist_step].complete_filename(self.file_loc)
            outputs_to_check.append(step_object_dict[worklist_step].full_filename)
        # $if should only fail if an extract only, but that may occur
        if outputs_to_check != []:
            proceed = self.check_file_removal(outputs_to_check)
            if not proceed:
                return
        
        # $now we need to get input and do the work. each step can only occur
        # $once and they occur in order. so we will write them in order
        # todo$ see if we can compress the code and make sure it is readable
        previous_output_file = ""
        extracted_files = []
        make_table_in_order = True
        for analysis_step in worklist:
            if analysis_step == "Extract":
                # $no if for this one, if extract is here it is the start
                id_file = self.collect_single_file("ID", "Extract", "CSV (*.csv)")
                if id_file == "": return
                # $always check if is good since it is first
                infile_is_good = self.check_input(step_object_dict[analysis_step],
                                                  id_file)
                if not infile_is_good:  return
                #$infile_is_good is just a check for blanks in the input file
                data_is_good = self.check_extractor_input(id_file,required_data_extractor_data, autofill_columns)
                if not data_is_good:  return
                mzml_files = self.collect_multiple_files("Centroided Data",
                                                         analysis_step, "mzML (*.mzML)")
                if mzml_files == []: return
                
                mzml_filenames = [os.path.basename(filename) for filename in
                                  mzml_files]
                extracted_files = [filename.replace(".mzML", ".tsv") for
                                   filename in mzml_filenames]
                extracted_files = [os.path.join(output_folder, filename) for
                                   filename in extracted_files]
                extracted_intermediate_files = extracted_files
                
                needed_files = list(set(extracted_files + extracted_intermediate_files))
                proceed = self.check_file_removal(needed_files)
                if not proceed:  return
                # $need to run the table if necessary. taken from the
                # $"Provide Time and Enrichment" elif
                if "Provide Time and Enrichment" in worklist:
                    previous_output_file = step_object_dict[
                        "Provide Time and Enrichment"].full_filename
                    self.get_data_table = TimeEnrichmentWindow(self,
                                                               extracted_files, previous_output_file)
                    self.get_data_table.exec_()
                    # $don't make the table twice
                    make_table_in_order = False
                    # $now that the table is done we need to confirm the user
                    # $hit the proceed button on the table (same check as in
                    # $elif analysis_step == "Theory Generation" )
                    if not os.path.exists(previous_output_file): return
                # $ modified from the extract-dir argument from the command line
                for m in tqdm(range(len(mzml_files)), total=len(mzml_files), desc="Extracting mzml files: "):
                    extractor = Extractor(
                        id_path=os.path.join(self.file_loc, id_file),
                        mzml_path=mzml_files[m],
                        out_path=extracted_intermediate_files[m],
                        settings_path=rate_settings_file,
                    )
                    extractor.load()
                    extractor.run()
                    extractor.write()
                    #$need to delete classes when they're done or they may linger in RAM
                    del extractor
            elif analysis_step == "Provide Time and Enrichment" and make_table_in_order:
                # $if coming right after a list
                if extracted_files == []:
                    extracted_files = self.collect_multiple_files(
                        "Extracted Data",
                        "Provide Time and Enrichment",
                        "TSV (*.tsv)"
                    )
                    if extracted_files == []: return
                    # $ensure the input files are good. only need to deal with
                    # $if the user just selected
                    for e_file in extracted_files:
                        infile_is_good = self.check_input(
                            step_object_dict[analysis_step],
                            e_file)
                        if not infile_is_good: return
                
                # $ now that we have the extracted files we can make a table
                # $the talbe will handle the output
                previous_output_file = step_object_dict[
                    analysis_step].full_filename
                self.get_data_table = TimeEnrichmentWindow(self,
                                                           extracted_files, previous_output_file)
                self.get_data_table.exec_()
            elif analysis_step == "Theory Generation":
                # $since the files are in the table can just read that in
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "time and enrichment",
                        analysis_step,
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step],
                        previous_output_file)
                    if not infile_is_good: return
                # $else is to deal with a failed write from the previous table
                # $ don't need an error message just return
                elif not os.path.exists(previous_output_file):
                    return
                
                # $final check to see if all of the files in the input table
                # $still exist.  don't want to error out in the middle of
                # $multiprocessing
                final_proceed = self.check_files_from_files(
                    previous_output_file, 0)
                if not final_proceed: return
                
                theorist = TheoryPreparer(
                    enrichment_path=previous_output_file,
                    out_path=step_object_dict[analysis_step].full_filename,
                    settings_path=rate_settings_file
                )
                theorist.prepare()
                theorist.write()
                del theorist
                previous_output_file = step_object_dict[analysis_step].full_filename
            elif analysis_step == "Fraction New Calculation":
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "theoretical output",
                        analysis_step,
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return
                    
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step],
                        previous_output_file)
                    if not infile_is_good: return
                # $not sure why this would happen but we'll put it here
                # $to avoid future error
                elif not os.path.exists(previous_output_file):
                    return
                fnewcalc = FractionNewCalculator(
                    model_path=previous_output_file,
                    out_path=step_object_dict[analysis_step].full_filename,
                    settings_path=rate_settings_file
                )
                fnewcalc.generate()
                if fnewcalc.error != "":
                    QtWidgets.QMessageBox.information(self, "Error", fnewcalc.error)
                    return
                fnewcalc.write()
                del fnewcalc
                previous_output_file = step_object_dict[
                    analysis_step].full_filename
            elif analysis_step == "Rate Calculation":
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "fraction new",
                        analysis_step,
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return
                    # $need to ensure that we have proper colmns which varies
                    # $by setting
                    needed_columns = [
                        settings.peptide_analyte_id_column,
                        settings.peptide_analyte_name_column,
                        "sample_group"]

                    if settings.use_abundance != "No":
                        needed_columns.extend(["abund_fn", "frac_new_abunds_std_dev"])
                    if settings.use_neutromer_spacing:
                        needed_columns.extend(["nsfn", "frac_new_mzs_std_dev"])
                    if settings.use_abundance != "No" and settings.use_neutromer_spacing:
                        needed_columns.extend(["cfn", "frac_new_combined_std_dev"])
                    step_object_dict[analysis_step].required_columns = needed_columns
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step],
                        previous_output_file)
                    if not infile_is_good: return
                
                # $need to get a graph folder and ensure it exists
                # $don't worry about overwriting files
                GraphFolder = os.path.join(self.file_loc, "Graph_Folder")
                MainGuiObject._make_folder(GraphFolder)
                ratecalc = RateCalculator(
                    model_path=previous_output_file,
                    out_path=step_object_dict[analysis_step].full_filename,
                    graph_folder=GraphFolder,
                    settings_path=rate_settings_file
                )
                ratecalc.calculate()
                ratecalc.write()
                del ratecalc
        QtWidgets.QMessageBox.information(self, "Success",
                                          "Analysis completed successfully")
    
    def check_table_checklist(self):
        current_worklist = []
        error_check = []
        for i in range(self.RequestedStepsTable.rowCount()):
            if self.RequestedStepsTable.item(i, 0).checkState() == QtCore.Qt.Checked:
                current_worklist.append(str(
                    self.RequestedStepsTable.item(i, 0).text()))
                error_check.append(i)
        # $the point of the error check is to ensure ther are no gaps in the
        # $checklist since then we will be missing critical info for the next
        # $step
        proper_length = max(error_check) - min(error_check) + 1
        if len(current_worklist) < proper_length:
            return ("There are gaps in your worklist.  Please check all boxes "
                    "between the first and last checked box")
        return current_worklist
    
    # $quick way to get a response to a question to the user.  exiting is a "No"
    def question_for_user(self, message):
        response = QtWidgets.QMessageBox.question(self, "Question", message,
                                                  QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if response == QtWidgets.QMessageBox.Yes:
            return True
        else:
            return False
    
    def large_text_question_for_use(self, title, infoText, detailedText):
        question = QtWidgets.QMessageBox(self)
        question.setWindowTitle("Question")
        question.setText(title)
        question.setIcon(QtWidgets.QMessageBox.Question)
        
        question.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        question.setDefaultButton(QtWidgets.QMessageBox.No)
        
        question.setDetailedText(detailedText)
        question.setInformativeText(infoText)
        
        question.setStyleSheet("QMessageBox{min-width:650 px;}")
        response = question.exec_()
        
        if response == QtWidgets.QMessageBox.Yes:
            return True
        else:
            return False
    
    def large_text_information(self, title, text, detailedText):
        info = QtWidgets.QMessageBox(self)
        info.setWindowTitle(title)
        info.setText(text)
        info.setIcon(QtWidgets.QMessageBox.Information)
        
        info.setStandardButtons(QtWidgets.QMessageBox.Ok)
        info.setDefaultButton(QtWidgets.QMessageBox.Ok)
        
        info.setDetailedText(detailedText)
        
        info.setStyleSheet("QMessageBox{min-width:650 px;}")
        info.exec_()
    
    #$ask user for one file
    def collect_single_file(self, to_load, step_name, load_type):
        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
                "file to load for {} step".format(to_load, step_name)))
        filename, file_type = QtWidgets.QFileDialog.getOpenFileName(self, 
                "Choose {} file to load".format(to_load), self.file_loc, 
                load_type)
        return filename
    
    #$ask user for potentially more than one file
    def collect_multiple_files(self, to_load, step_name, load_type):
        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
                "files to load for {} step".format(to_load, step_name)))
        filenames, file_type = QtWidgets.QFileDialog.getOpenFileNames(self, 
                "Choose {} file to load".format(to_load), self.file_loc, 
                load_type)
        return filenames
    
    def change_settings(self):
        self.set_menu = rate_settings.Rate_Setting_Menu(self, rate_settings_file)
        self.set_menu.show()
    
    def change_converter_settings(self):
        self.set_menu = guide_settings.Converter_Setting_Menu(self, id_settings_file)
        self.set_menu.show()
    
    #$checks if a file exists, and if we can actually do so.  if we can't
    #$an error will be triggered later so we'll sort that out here.  
    #$has the advantage of deleting the files so the user can't open a file between now and when we write it
    #$ask the user about this for everything but settings (that is dealt with better in the run_rate_worklist function)
    def check_file_removal(self, list_of_filenames, ask_permission = True):
        files_to_remove = []
        open_files = []
        #$find files that exist
        for filename in list_of_filenames:
            if os.path.exists(filename):
                files_to_remove.append(filename)
        # $let the user decide if we should continue.
        if files_to_remove != []:
            if ask_permission:
                proceed = self.large_text_question_for_use("Some files already exist and will be overwritten.",
                                                           "Do you still wish to proceed?",
                                                           "Files to be overwritten:\n" + ",\n".join(files_to_remove))
                if not proceed:
                    return False
            for filename in files_to_remove:
                try:
                    os.remove(filename)
                except PermissionError:
                    open_files.append(filename)
            if open_files != []:
                self.large_text_information("Error", "Some files cannot be overwritten.\n\n "
                                                     "They are likely open in another program. Please close "
                                                     "and try again.",
                                            "Files unable to be opened:\n" + ",\n".join(open_files))
    
                return False
        # $will return true if no files already exist or the user wants to
        # $overwrite and they can be removed so we have permission
        return True
    
    def check_input(self, relevant_object, filename):
        has_needed_columns = relevant_object.check_input_file(filename)
        if not has_needed_columns:
            QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                                                              "missing needed columns. Please correct and try again".format(
                filename)))
        return has_needed_columns
    
    # $there are cases (specifically the table going to theory) where there
    #$is a chance that the files referenced in a id file may not exist
    # $this will check for that
    def check_files_from_files(self, input_file, filename_column):
        with open(input_file, 'r') as infile:
            if input_file[-4:] == ".tsv":
                reader = csv.reader(infile, delimiter="\t")
            elif input_file[-4:] == ".csv":
                reader = csv.reader(infile)
            next(reader)  # $skip header
            # $may wish to add a check that rows are of appropriate length
            for row in reader:
                if not os.path.exists(row[filename_column]):
                    QtWidgets.QMessageBox.information(self, "Error", ("File {} "
                                                                      "could not be found. Please correct input file and try "
                                                                      "again".format(row[filename_column]
                                                                                     )))
                    return False
        return True
    
    # $ensures a folder exists by making it if it does not.
    @staticmethod
    def _make_folder(folder):
        if not os.path.isdir(folder):
            os.makedirs(folder)
            """
            #$this is for is we want to delete the entire folder.  this particular use is far too crude.  Any subfolders or files
            #$even unrelated to DeuteRater will be deleted.  If we need to put this back in we need to target only DeuteRater files
            #$that should be unnecessary because we have the the overwrite warning
            if os.path.isdir(folder):
                answer = self.question_for_user(f"Would you like to delete the contents of folder {folder}?")
                if answer:
                    rmtree(folder)
                    os.makedirs(folder)
            else:
                os.makedirs(folder)"""
        else:
            if os.path.isdir(folder):
                rmtree(folder)
            os.makedirs(folder)
    
    @staticmethod
    def _get_file_names(folder_loc, operations, object_dict):
        needed_files = []
        for o in operations:
            if o != "Extract":
                full_filename = os.path.join(folder_loc,
                                             object_dict.output_filename)
                needed_files.append(full_filename)
        return needed_files

    #$since the users will mainly be filling in a template we need to check the input.  we need to check that the appropriate columns have data in them
    #$if there is something that we can autofill, we'll do that in the extractor itself since that is already reading in the file
    def check_extractor_input(self, filename, needed_columns, autofill_columns):
        df = pd.read_csv(filename)
        error_message = "Column \"{}\" has blank cells within it.  Please correct and try again."
        for n in needed_columns:
            #$ blanks will be true, so this is number of blanks
            #$works for text columns
            if sum(df[n] == "") > 0:
                QtWidgets.QMessageBox.information(self, "Error", error_message.format(n))
                return False
            #$numerical columns have nan's not blanks
            elif sum(df[n].isnull()) > 0:
                QtWidgets.QMessageBox.information(self, "Error", error_message.format(n))
                return False
        
        for a in autofill_columns:
            if sum(df[a].isnull()) > 0 or sum(df[a] == "") > 0:
                perform_autofill = True
                bad_column = a
                break
        else:
            perform_autofill = False
        if perform_autofill:
            result = self.autofill(df, filename, bad_column)
            if not result:
                return result
            
        return True
    
    #$this is based on _interpret_aa_sequences from peaksXPro
    #$if the user has not provided a chemical formula, literature n value and theoretical neutral mass
    #$for now we'll just assume that if any are blank we should just autofill everything
    def autofill(self, df, guide_file_location, bad_column):
        converter_settings.load(id_settings_file)
        aa_comp_df = pd.read_csv(converter_settings.aa_elem_comp_path, sep='\t')
        aa_comp_df.set_index('amino_acid', inplace=True)

        aa_label_df = pd.read_csv(converter_settings.aa_label_path, sep='\t')
        aa_label_df.set_index('study_type', inplace=True)
        # TODO: Find out where to store settings, then decide which 'studytype'
        #       to use as default.
        aa_labeling_dict = aa_label_df.loc[converter_settings.study_type, ].to_dict()

        elem_df = pd.read_csv(converter_settings.elems_path, sep='\t')
        # This is necessary if we have all of the different isotopes in the tsv
        element_index_mask = [
            0,  # Hydrogen
            10,  # Carbon-12
            13,  # Nitrogen-14
            15,  # Oxygen-16,
            30, # Phosphorous
            31  # Sulfer-32
        ]
        elem_df = elem_df.iloc[element_index_mask]
        elem_df.set_index('isotope_letter', inplace=True)
        #$if dtypes are wrong we can either error out or have improper values, so we'll force it
        df['cf'] = df['cf'].astype(str) 
        df['Peptide Theoretical Mass'] = df['Peptide Theoretical Mass'].astype(float)
        df['literature_n'] = df['literature_n'].astype(float)
        for row in df.itertuples():
            i = row.Index
            aa_counts = {}
            for aa in row.Sequence:
                if aa not in aa_counts.keys():
                    aa_counts[aa] = 0
                aa_counts[aa] += 1
            elem_dict = peputils.calc_cf(aa_counts, aa_comp_df)
            theoretical_mass = peputils.calc_theory_mass(elem_dict, elem_df)
            literature_n = peputils.calc_add_n(aa_counts, aa_labeling_dict)
            df.at[i, 'cf'] = ''.join(
                k + str(v) for k, v in elem_dict.items() if v > 0
            )
            df.at[i, 'Peptide Theoretical Mass'] = theoretical_mass
            df.at[i, 'literature_n'] = literature_n
        #$now we need to srite out.  since we only needed read permissions elsewhere we may have a problem.
        #$but it does ensure we need to pass a df around.
        try:
            df.to_csv(guide_file_location, index = False)
            return True
        except IOError:
                QtWidgets.QMessageBox.information(self, "Error", 
                    ("Column \"{}\" in file \"{}\" needed to be filled in by DeuteRater-H."
                     " This file is open in another program or DeuteRater-H does not have "
                     "write permission to this location. Please either fill in the column or close "
                     "the file and try again.".format(bad_column, guide_file_location)))
                return False

# $since we have to load settings in each file, and need a way to adjust
# $settings, we'll
def make_temp_file(filename, new_filename):
    copyfile(filename, new_filename)

def main():
    # $needed for windows multiprocessing which will happen at some point
    import sys
    
    mp.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    gui_object = MainGuiObject(None)
    gui_object.show()
    app.exec_()

if __name__ == '__main__':
    main()
