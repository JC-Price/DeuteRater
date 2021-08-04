# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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


#$similar to deuterater step the goal is to force the the correct headers.
#$useful for peaks only. could adjust this and deuterater step to inherit from a common parent class
#$or merge into one class, but not worth it right now
class deuteconvert_peaks_required_headers(object):
    def __init__(self,protein_header, protein_peptide_header, features_header):
        self.protein_header = protein_header
        self.protein_peptide_header = protein_peptide_header
        self.features_header = features_header
        
    def protein_file_check(self, full_path):
        return self.check_input_file(full_path, self.protein_header)
    
    def protein_peptide_check(self, full_path):
        return self.check_input_file(full_path, self.protein_peptide_header)
    
    def features_check(self, full_path):
        return self.check_input_file(full_path, self.features_header)
        
    def check_input_file(self, full_path, header):
        with open (full_path, 'r') as infile:
            if full_path[-4:] == ".tsv":
                reader = csv.reader(infile, delimiter = "\t")
            elif full_path[-4:] == ".csv":
                reader = csv.reader(infile)
            else:
                raise ValueError("Improper input file")
            firstrow = next(reader)
        return set(header).issubset(firstrow)

#$this class holds data for each deuterater step
#$allows easier calls in the main and one large batch of declarations
#$limits hard coding issues
class deuterater_step(object):
    def __init__(self, output_filename, peptide_required_columns,
                 lipid_required_columns):
        self.output_filename = output_filename
        self.peptide_required_columns = peptide_required_columns
        self.lipid_required_columns = lipid_required_columns
        
    def complete_filename(self, folder):
        self.full_filename = os.path.join(folder, self.output_filename)
        
    def check_input_file(self, input_file, biomolecule_type):
        with open (input_file, 'r') as infile:
            if input_file[-4:] == ".tsv":
                reader = csv.reader(infile, delimiter = "\t")
            elif input_file[-4:] == ".csv":
                reader = csv.reader(infile)
            else:
                raise ValueError("Improper input file")
            firstrow = next(reader)
        if biomolecule_type == "Peptide":
            return set(self.peptide_required_columns).issubset(firstrow)
        elif biomolecule_type == "Lipid":
            return set(self.lipid_required_columns).issubset(firstrow)
    
#$ we're going to make two classes. these will be here so we can swap between
#$two different settings menus as needed.  we need string and numerical classes
#$technically we could have a parent and two children, but we'll see if we have
#$enough similarities to be worth it
class setting_numerical_info(object):
    #$data_object is the object from the menu, starting value is the initial 
    #$value, setting_name is tha name in the settings menu, integer is a 
    #$boolean to determine if it is an integer (True) or a float (False)
    def __init__(self, data_object, setting_name, starting_value, integer):
        self.data_object = data_object
        self.current_value = starting_value
        self.setting_name = setting_name
        self.integer = integer
        
    def set_object_value(self):
        self.data_object.setValue(self.current_value)
    
    def save_value(self):
        if self.integer:
            self.current_value = int(self.data_object.value())
        else:
            self.current_value = float(self.data_object.value())
        return self.setting_name, self.current_value
    #$the purpose of this is just to compare if the values match without
    #$saving so we can compare things when closing
    def compare_value(self):
        if self.integer:
            return self.current_value == int(self.data_object.value())
        else:
            return self.current_value == float(self.data_object.value())

class setting_string_info(object):
    #$arguments are the same as setting_numerical_info except true_false
    #$true_false determines if it is a yes no(True) or just a string
    def __init__(self, data_object, setting_name, starting_value, true_false):
        self.data_object = data_object
        self.current_value = starting_value
        self.setting_name = setting_name
        self.true_false = true_false
        
    def set_object_value(self):
        if self.true_false:
            if self.current_value:
                string_value = "Yes"
            else:
                string_value = "No"
        else:
            string_value = self.current_value
        index = self.data_object.findText(string_value)
        self.data_object.setCurrentIndex(index)
        
    def save_value(self):
        value = str(self.data_object.currentText())
        if self.true_false:
            if value == "Yes":
                self.current_value = True
            else:
                self.current_value = False
        else:
            self.current_value = value
        return self.setting_name, self.current_value
    
    #$the purpose of this is just to compare if the values match without
    #$saving so we can compare things when closing
    def compare_value(self):
        value = str(self.data_object.currentText())
        if self.true_false:
            if value == "Yes" and self.current_value:
                return True
            elif value == "No" and not self.current_value:
                return True
            else:
                return False
        else:
            return value == self.current_value

class setting_checkbox_info(object):
    def __init__(self, data_object, setting_name, starting_value):
        self.data_object = data_object
        self.current_value = starting_value
        self.setting_name = setting_name
    
    def set_object_value(self):
        self.data_object.setChecked(self.current_value)
    
    def save_value(self):
        self.current_value = self.data_object.isChecked()
        return self.setting_name, self.current_value
    
    def compare_value(self):
        return self.current_value == self.data_object.isChecked()
