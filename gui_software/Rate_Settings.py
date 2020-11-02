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

import deuterater.settings as settings
from utils.useful_classes import setting_numerical_info, setting_string_info


#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

settings_file_ui_location = os.path.join(location, "ui_files", "Settings_Menu.ui")
loaded_ui = uic.loadUiType(settings_file_ui_location)[0]

class Rate_Setting_Menu(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, current_setting_file = None):
        super(Rate_Setting_Menu, self).__init__(parent)
        settings.load(current_setting_file)
        self.current_setting_file =current_setting_file
        #$this is needed to slim things down a bit
        self.setWindowTitle("Rate Settings Menu")
        self.setupUi(self)
        self.all_settings=[
            setting_string_info(self.recognize_available_cores, "recognize_available_cores",
                                 settings.recognize_available_cores, True),
            setting_numerical_info(self.default_cores, "n_processors",
                                   settings.n_processors, True),
            setting_string_info(self.rt_unit, "id_file_rt_unit",
                                settings.id_file_rt_unit, False),
            setting_numerical_info(self.time_window, "time_window",
                                    settings.time_window, False),
            setting_numerical_info(self.ppm_error, "ppm_window",
                                    settings.ppm_window, True),
            setting_string_info(self.heavy_label, "heavy_isotope",
                                settings.heavy_isotope, False),
            setting_string_info(self.calculate_n_values, "use_empir_n_value",
                                settings.use_empir_n_value, True),
            setting_string_info(self.use_abundance, "use_abundance",
                                settings.use_abundance, True),
            setting_string_info(self.use_neutromer_spacing, "use_neutromer_spacing",
                                settings.use_neutromer_spacing, True),
            setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                   settings.minimum_nonzero_points, True),
            setting_string_info(self.roll_up_option,"roll_up_rate_calc",
                                 settings.roll_up_rate_calc, True),
            setting_string_info(self.asymptope_type, "asymptote", 
                                 settings.asymptote, False),
            setting_numerical_info(self.fixed_asymptote_value, "fixed_asymptote_value",
                                   settings.fixed_asymptote_value, False),
            setting_numerical_info(self.proliferation_adjustment, "proliferation_adjustment",
                                   settings.proliferation_adjustment, False),
            setting_string_info(self.bias_selection_option, "bias_calculation",
                                 settings.bias_calculation, False),
            setting_numerical_info(self.abund_manual_bias, "abundance_manual_bias",
                                   settings.abundance_manual_bias, False),
            setting_numerical_info(self.spacing_manual_bias, "spacing_manual_bias",
                                   settings.spacing_manual_bias, False),
            setting_numerical_info(self.combined_manual_bias, "combined_manual_bias",
                                   settings.combined_manual_bias, False),
            setting_numerical_info(self.min_allowed_m0_change, "min_allowed_abund_max_delta",
                                   settings.min_allowed_abund_max_delta, False),
            setting_numerical_info(self.min_sequence_length, "min_aa_sequence_length",
                                   settings.min_aa_sequence_length, True),
            setting_numerical_info(self.min_n_value, "min_allowed_n_values", 
                                   settings.min_allowed_n_values, True)
            ]
        for setting_object in self.all_settings:
            setting_object.set_object_value()
        self.SaveButton.clicked.connect(self.save_settings)
        self.ExitButton.clicked.connect(self.close)
        
    def save_settings(self):
        #$we need to provide the values that are not altred for the dump
        save_value_dict = Rate_Setting_Menu._get_filters()
        for setting_object in self.all_settings:
            name, value = setting_object.save_value()
            save_value_dict[name] = value
        settings.freeze(self.current_setting_file, save_value_dict)
        
    def check_for_changes(self):
        for setting in self.all_settings:
            if not setting.compare_value():
                return False
        return True
            
    
    #$should overwrite the close of the exit button and the red x in the corner  
    def closeEvent(self, event):
        if self.check_for_changes():
            event.accept()
        else:
            reply = QtWidgets.QMessageBox.question(self, "Unsaved Changes", 
                ("There are unsaved changes.  Would you like to save before"
                " exiting?"), 
                QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No|QtWidgets.QMessageBox.Cancel)
            if reply == QtWidgets.QMessageBox.No:
                event.accept()        
            elif reply == QtWidgets.QMessageBox.Yes:
                self.save_settings()
                event.accept()
            #$hit red x or cancel just don't exit
            else: 
                event.ignore()  
      
    #$the point of this is to get the values from the settings .yaml 
    #$that the user is not altering. don't need special classes since 
    #$ we have to officially declare and all we need to do is get the value
    @staticmethod
    def _get_filters():
        unalterable_settings = {
            "debug_level" : settings.debug_level,
            "trim_ids_to_mzml_bounds" : settings.trim_ids_to_mzml_bounds,
            "fpe_tolerance" : settings.fpe_tolerance,
            "chunk_size" : settings.chunk_size,
            "chunking_method_threshold" : settings.chunking_method_threshold,
            "max_valid_angle" : settings.max_valid_angle,
            "peak_ratio_denominator" : settings.peak_ratio_denominator,
            "peptide_analyte_id_column" : settings.peptide_analyte_id_column,
            "lipid_analyte_id_column" : settings.lipid_analyte_id_column,
            "peptide_analyte_name_column" : settings.peptide_analyte_name_column,
            "lipid_analyte_name_column" : settings.lipid_analyte_name_column,
            "unique_sequence_column" : settings.unique_sequence_column,
            "maximum_theoretical_pct" : settings.maximum_theoretical_pct,
            "labeling_step_size" : settings.labeling_step_size,
            "peak_lookback" : settings.peak_lookback,
            "peak_lookahead" : settings.peak_lookahead,
            "baseline_lookback" : settings.baseline_lookback,
            "min_envelopes_to_combine" : settings.min_envelopes_to_combine,
            "zscore_cutoff": settings.zscore_cutoff,
            "mz_proximity_tolerance" : settings.mz_proximity_tolerance,
            "error_of_zero" : settings.error_of_zero,
            "abundance_agreement_filter" : settings.abundance_agreement_filter,
            "spacing_agreement_filter" : settings.spacing_agreement_filter,
            "combined_agreement_filter" : settings.combined_agreement_filter,
            "error_of_non_replicated_point" : settings.error_of_non_replicated_point,
            "y_intercept_of_fit" : settings.y_intercept_of_fit,
            "enrichement_of_zero" : settings.enrichement_of_zero,
            "minimum_abund_change": settings.minimum_abund_change
            }
        return unalterable_settings
        