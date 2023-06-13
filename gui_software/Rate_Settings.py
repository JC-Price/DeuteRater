# -*- coding: utf-8 -*-
"""
Copyright (c) 2016-2020 Bradley Naylor,  J.C. Price, and Brigham Young University
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
import pandas as pd
from PyQt5 import uic, QtWidgets

import rater.settings as settings
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
        
        self.fill_study_type_combobox()
        self.all_settings=[
            setting_string_info(self.recognize_available_cores, "recognize_available_cores",
                                settings.recognize_available_cores, True),
            setting_numerical_info(self.default_cores, "n_processors",
                                   settings.n_processors, True),
            setting_string_info(self.study_type_combobox, "study_type",
                                settings.study_type, False),
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
                                settings.use_abundance, False),
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
                                   settings.min_allowed_n_values, True),
            setting_numerical_info(self.ms_level, "ms_level",
                                   settings.ms_level, True),
            setting_string_info(self.use_chromatography_division, "use_chromatography_division",
                                settings.use_chromatography_division, False),
            setting_string_info(self.verbose_rate, "verbose_rate",
                                settings.verbose_rate, True)
            ]
        for setting_object in self.all_settings:
            setting_object.set_object_value()
        self.SaveButton.clicked.connect(self.save_settings)
        self.ExitButton.clicked.connect(self.close)
        self.LoadButton.clicked.connect(self.load_settings)
        self.setWindowTitle("Rate Calculator Settings")
        
    def fill_study_type_combobox(self):
        temp_df = pd.read_csv(settings.aa_label_path, sep ="\t")
        self.study_types =  list(temp_df["study_type"].unique())
        for study_type in self.study_types:
            self.study_type_combobox.addItem(study_type)
        
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
        
    def load_settings(self):
        response = QtWidgets.QMessageBox.question(self, "Question", "Would you like to load a already existing settings file? This will overwrite all current settings.",
                                                  QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if response == QtWidgets.QMessageBox.No:
            return

        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the settings "
                                                         "file to load"))
        filename, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
                                                                    "Choose settings file to load",
                                                                    os.path.curdir,
                                                                    "*.yaml",
                                                                    options=QtWidgets.QFileDialog.DontUseNativeDialog)
        
        comp_results = settings.compare(self.current_setting_file, filename)
        if comp_results == "Error":
            QtWidgets.QMessageBox.warning(self, "Error", ("Issue reading .yaml file. Please make sure the .yaml still exists and is not currently opened."))
            return
        elif comp_results == "Different Keys":
            QtWidgets.QMessageBox.warning(self, "Error", ("Loaded settings file either is missing settings or has too many. Please try a different file with the correct settings"))
            return
        elif comp_results == "MATCH":
            QtWidgets.QMessageBox.information(self, "Info", ("Settings file have the same data."))
            return
        elif comp_results == "Mismatched Keys":
            settings.load(filename)
            settings.freeze(self.current_setting_file)
            self.all_settings = [
                setting_string_info(self.recognize_available_cores, "recognize_available_cores",
                                    settings.recognize_available_cores, True),
                setting_numerical_info(self.default_cores, "n_processors",
                                       settings.n_processors, True),
                setting_string_info(self.study_type_combobox, "study_type",
                                    settings.study_type, False),
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
                                    settings.use_abundance, False),
                setting_string_info(self.use_neutromer_spacing, "use_neutromer_spacing",
                                    settings.use_neutromer_spacing, True),
                setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                       settings.minimum_nonzero_points, True),
                setting_string_info(self.roll_up_option, "roll_up_rate_calc",
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
                                       settings.min_allowed_n_values, True),
                setting_numerical_info(self.ms_level, "ms_level",
                                       settings.ms_level, True),
                setting_string_info(self.use_chromatography_division, "use_chromatography_division",
                                    settings.use_chromatography_division, False),
                setting_string_info(self.verbose_rate, "verbose_rate",
                                    settings.verbose_rate, True)
            ]
            for setting_object in self.all_settings:
                setting_object.set_object_value()
            QtWidgets.QMessageBox.information(self, "Info", ("Settings successfully loaded."))
            return

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
            'study_type': settings.study_type,
            "peak_ratio_denominator" : settings.peak_ratio_denominator,
            "peptide_analyte_id_column" : settings.peptide_analyte_id_column,
            "lipid_analyte_id_column" : settings.lipid_analyte_id_column,
            "peptide_analyte_name_column" : settings.peptide_analyte_name_column,
            "aa_labeling_sites_path": settings.aa_label_path,
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
            "minimum_abund_change": settings.minimum_abund_change,
            "intensity_filter": settings.intensity_filter,
            "rel_height": settings.rel_height,
            "sampling_rate": settings.sampling_rate,
            "smoothing_width": settings.smoothing_width,
            "smoothing_order": settings.smoothing_order,
            "allowed_peak_variance_min": settings.allowed_peak_variance_min,
            "adduct_weight": settings.adduct_weight,
            "variance_weight": settings.variance_weight,
            "ID_weight": settings.ID_weight,
            "intensity_weight": settings.intensity_weight,
            "how_divided": settings.how_divided,
            "allowed_neutromer_peak_variance": settings.allowed_neutromer_peak_variance,
            "rate_output_format": settings.rate_output_format,
            "s_n_filter": settings.s_n_filter,
            'remove_filters': settings.remove_filters,
            "separate_adducts": settings.separate_adducts,
            "vertical_gridlines": settings.vertical_gridlines,
            }
        return unalterable_settings
