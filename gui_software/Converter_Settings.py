# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 08:04:31 2020

@author: Naylor

for the moment we are going to create two settings menus.  this is for the
sake of quick creation and separating the values
Will not use all of the settings, can swap in or out as necessary
Also since .yaml files can be opened by notepad and are not encrypted
we will not give the user the option to change the defaults in the gui.
we'll start with the rate settings as I am more familiar with
"""

"""
the following are settings that are in use in DeuteRater but are not alerable
in this menu.  if we need them just grab them
debug_level - adjusting debug setting should not be set in a gui
trim_ids_to_mzml_bounds - should not be needed often enough to put here
fpe_tolerance - should not be needed often enough to put here
chunk_size - should not be needed often enough to put here
chunking_method_threshold - should not be needed often enough to put here
max_valid_angle - should not be needed often enough to put here
peak_ratio_denominator - should not be needed often enough to put here
analyte_id_column - don't need to adjust until we have different id files
analyte_name_column - don't need to adjust until we have different id files
unique_sequence_column - don't need to adjust until we have different id files
maximum_theoretical_pct - not using currently
labeling_step_size - not using currently
peak_lookback - shouldn't need to alter often
peak_lookahead - shouldn't need to alter often
baseline_lookback - shouldn't need to alter often
min_envelopes_to_combine  - shouldn't need to alter often
zscore_cutoff- shouldn't need to alter often.
mz_proximity_tolerance- shouldn't need to alter often
error_of_zero - shouldn't need to alter often
abundance_agreement_filter - shouldn't need to adjust all that often
spacing_agreement_filter - shouldn't need to adjust all that often
combined_agreement_filter - shouldn't need to adjust all that often
error_of_non_replicated_point - roll up is rare enough that we should
    not need often
y_intercept_of_fit - do need to be adjustable but should not be usual
enrichement_of_zero - adjusting this is usually unnecessary. troubleshooting
        only so no need for normal settings

"""
import os

from PyQt5 import uic, QtWidgets

import deuteconvert.settings as settings
from utils.useful_classes import setting_numerical_info, setting_string_info, setting_checkbox_info

# location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

settings_file_ui_location = os.path.join(location, "ui_files", "Converter_Settings_Menu.ui")
loaded_ui = uic.loadUiType(settings_file_ui_location)[0]


class Converter_Setting_Menu(QtWidgets.QDialog, loaded_ui):
	def __init__(self, parent=None, current_setting_file=None):
		super(Converter_Setting_Menu, self).__init__(parent)
		settings.load(current_setting_file)
		self.current_setting_file = current_setting_file
		# $this is needed to slim things down a bit
		self.setWindowTitle("Converter Settings Menu")
		self.setupUi(self)
		self.all_settings = [
			setting_numerical_info(self.min_charge, "min_charge_state", settings.min_charge_state, True),
			setting_numerical_info(self.max_charge, "max_charge_state", settings.max_charge_state, True),
		]
		for setting_object in self.all_settings:
			setting_object.set_object_value()
		self.SaveButton.clicked.connect(self.save_settings)
		self.ExitButton.clicked.connect(self.close)
	
	def save_settings(self):
		# Check if the minimum and maximum charge numbers makes since:
		if int(self.min_charge.value()) > int(self.max_charge.value()):
			QtWidgets.QMessageBox.information(self, "Error", ("The Minimum Charge must be less than or equal to "
													"the Maximum Charge. Settings are NOT saved. \n\nPlease make "
													"Minimum Charge less than or equal to Maximum Charge to save."))
			return False
		# $we need to provide the values that are not altred for the dump
		save_value_dict = Converter_Setting_Menu._get_filters()
		for setting_object in self.all_settings:
			name, value = setting_object.save_value()
			save_value_dict[name] = value
		settings.freeze(self.current_setting_file, save_value_dict)
		return True
	
	def check_for_changes(self):
		for setting in self.all_settings:
			if not setting.compare_value():
				return False
		return True
	
	# $should overwrite the close of the exit button and the red x in the corner
	def closeEvent(self, event):
		if self.check_for_changes():
			event.accept()
		else:
			reply = QtWidgets.QMessageBox.question(self, "Unsaved Changes",
												   ("There are unsaved changes.  Would you like to save before"
													" exiting?"),
												   QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)
			if reply == QtWidgets.QMessageBox.No:
				event.accept()
			elif reply == QtWidgets.QMessageBox.Yes:
				if self.save_settings():
					event.accept()
				event.ignore()
			# $hit red x or cancel just don't exit
			else:
				event.ignore()
			
			# $the point of this is to get the values from the settings .yaml
	
	# $that the user is not altering. don't need special classes since
	# $ we have to officially declare and all we need to do is get the value
	@staticmethod
	def _get_filters():
		unalterable_settings = {
			'mass_cutoffs': settings.mass_cutoffs,
			'rt_proximity_tolerance': settings.rt_proximity_tolerance,
			'mz_proximity_tolerance': settings.mz_proximity_tolerance,
			'start_time': settings.start_time,
			'study_type': settings.study_type,
			'aa_elemental_composition_path': settings.aa_elem_comp_path,
			'aa_labeling_sites_path': settings.aa_label_path,
			'elements_path': settings.elems_path,
			'post_translational_modifications_path': settings.ptms_path,
			'remove_duplicates': settings.remove_duplicates,
		}
		return unalterable_settings
