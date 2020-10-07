# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:52:01 2020

hold various classes that are just holding variables for easy reference
maybe one or two basic functions

@author: JCPrice
"""
import os
import csv

#$this class holds data for each deuterater step
#$allows easier calls in the main and one large batch of declarations
#$limits hard coding issues
class deuterater_step(object):
    def __init__(self, output_filename, required_columns):
        self.output_filename = output_filename
        self.required_columns = required_columns
        
    def complete_filename(self, folder):
        self.full_filename = os.path.join(folder, self.output_filename)
        
    def check_input_file(self, input_file):
        with open (input_file, 'r') as infile:
            if input_file[-4:] == ".tsv":
                reader = csv.reader(infile, delimiter = "\t")
            elif input_file[-4:] == ".csv":
                reader = csv.reader(infile)
            else:
                raise ValueError("Improper input file")
            firstrow = next(reader)
        return set(self.required_columns).issubset(firstrow)
  
    
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