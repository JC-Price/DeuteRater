# -*- coding: utf-8 -*-
"""
Copyright (c) 2024 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, Benjamin Driggs,
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

governs the settings menu including:
    loading current settings
    saving settings
    ensuring settings are valid
"""

# Todo: refactor load settings functions so we don't have duplicated code - Ben D

import os
import sys

from PyQt5 import uic, QtWidgets

import deuterater.settings as settings
from utils.useful_classes import setting_numerical_info, setting_string_info

# when compiling/building for an executable, set all of these to True, otherwise leave as False
# copy "exe_mode = False" and search using ctrl+shift+f to find each instance
exe_mode = False
if exe_mode:
    location = os.path.dirname(os.path.abspath(sys.executable))
else:
    location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

default_peptide_rate_settings = os.path.join(location, "resources", "peptide_settings.yaml")
default_lipid_rate_settings = os.path.join(location, "resources", "lipid_settings.yaml")

settings_file_ui_location = os.path.join(location, "ui_files", "Settings_Menu.ui")
loaded_ui = uic.loadUiType(settings_file_ui_location)[0]


class Rate_Setting_Menu(QtWidgets.QDialog, loaded_ui):

    def __init__(self, parent=None, current_setting_file=None):
        super(Rate_Setting_Menu, self).__init__(parent)
        settings.load(current_setting_file)
        self.current_setting_file = current_setting_file
        # this is needed to slim things down a bit
        self.setWindowTitle("Rate Settings Menu")
        self.setupUi(self)
        # initial set up of the different settings.
        self.all_settings = [
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
            setting_string_info(self.use_chromatography_division,
                                "use_chromatography_division",
                                settings.use_chromatography_division,
                                False),
            setting_string_info(self.label_key, "label_key",
                                settings.label_key, False),
            setting_numerical_info(self.min_AA_length, "min_aa_sequence_length",
                                   settings.min_aa_sequence_length, True),
            setting_numerical_info(self.min_allowed_n_value, "min_allowed_n_values",
                                   settings.min_allowed_n_values, True),
            setting_numerical_info(self.min_allowed_rate, "minimum_allowed_sequence_rate",
                                   settings.minimum_allowed_sequence_rate, False),
            setting_numerical_info(self.max_allowed_rate, "maximum_allowed_sequence_rate",
                                   settings.maximum_allowed_sequence_rate, False),
            setting_string_info(self.graph_file_type,
                                "graph_output_format",
                                settings.graph_output_format,
                                False),
            setting_numerical_info(self.abund_manual_bias, "abundance_manual_bias",
                                   settings.abundance_manual_bias, False),
            setting_string_info(self.asymptope_type, "asymptote",
                                settings.asymptote, False),
            setting_numerical_info(self.combined_manual_bias, "combined_manual_bias",
                                   settings.combined_manual_bias, False),
            setting_numerical_info(self.fixed_asymptote_value, "fixed_asymptote_value",
                                   settings.fixed_asymptote_value, False),
            setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                   settings.minimum_nonzero_points, True),
            setting_numerical_info(self.proliferation_adjustment, "proliferation_adjustment",
                                   settings.proliferation_adjustment, False),
            setting_string_info(self.use_outlier_removal, "use_outlier_removal",
                                settings.use_outlier_removal, True),
            setting_numerical_info(self.spacing_manual_bias, "spacing_manual_bias",
                                   settings.spacing_manual_bias, False),
            setting_string_info(self.abundance_type, "abundance_type",
                                settings.abundance_type, False),
            setting_string_info(self.calculate_n_values, "use_empir_n_value",
                                settings.use_empir_n_value, True),
            setting_string_info(self.verbose_rate, "verbose_rate",
                                settings.verbose_rate, True),
            setting_string_info(self.bias_calculation, "bias_calculation",
                                settings.bias_calculation, False),
            setting_numerical_info(self.max_fn_standard_deviation, "max_fn_standard_deviation",
                                   settings.max_fn_standard_deviation, False),
            setting_string_info(self.fraction_new_calculation, "fraction_new_calculation",
                                settings.fraction_new_calculation, False),
            setting_numerical_info(self.n_value_cv_limit, "n_value_cv_limit",
                                   settings.n_value_cv_limit, False),
            setting_string_info(self.separate_adducts, "separate_adducts",
                                settings.separate_adducts, True),
            setting_numerical_info(self.r2_threshold, "r2_threshold",
                                   settings.r2_threshold, False),
            setting_numerical_info(self.zscore_cutoff, "zscore_cutoff",
                                   settings.zscore_cutoff, False),
            setting_string_info(self.graph_n_value_calculations, "graph_n_value_calculations",
                                settings.graph_n_value_calculations, True),
            setting_string_info(self.save_n_value_data, "save_n_value_data",
                                settings.save_n_value_data, True)
        ]
        self.setWindowTitle("Settings Menu")
        for setting_object in self.all_settings:
            setting_object.set_object_value()

        self.LoadPeptideButton.clicked.connect(self.load_peptide_default_settings)
        self.LoadLipidButton.clicked.connect(self.load_lipid_default_settings)
        self.LoadButton.clicked.connect(self.load_settings)
        self.SaveButton.clicked.connect(self.save_settings)
        self.ExitButton.clicked.connect(self.close)

    # save the current settings
    def save_settings(self):
        # we need to provide the values that are not altered for the dump
        save_value_dict = Rate_Setting_Menu._get_filters()
        for setting_object in self.all_settings:
            name, value = setting_object.save_value()
            save_value_dict[name] = value

        settings.freeze(self.current_setting_file, save_value_dict)
        return True

    # if the user wants to exit we can check if something has changed, so we can give them a chance to save
    def check_for_changes(self):
        for setting in self.all_settings:
            if not setting.compare_value():
                return False
        return True

    # load a previous saved settings file
    def load_settings(self):
        response = QtWidgets.QMessageBox.question(self, "Question",
                                                  "Would you like to load a already existing settings file? This will overwrite all current settings.",
                                                  QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if response == QtWidgets.QMessageBox.No:
            return

        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the settings "
                                                         "file to load"))
        filename, file_type = QtWidgets.QFileDialog.getOpenFileName(self,
                                                                    "Choose settings file to load",
                                                                    location,
                                                                    "*.yaml",
                                                                    options=QtWidgets.QFileDialog.DontUseNativeDialog)
        # only bother with an error message if there is an error.  if the user just exited without selecting anything don't bother them with it
        if filename == "":
            return

        comp_results = settings.compare(self.current_setting_file, filename)
        if comp_results == "Error":
            QtWidgets.QMessageBox.warning(self, "Error", (
                "Issue reading .yaml file. Please make sure the .yaml still exists and is not currently opened."))
            return
        elif comp_results == "Different Keys":
            QtWidgets.QMessageBox.warning(self, "Error", (
                "Loaded settings file either is missing settings or has too many. Please try a different file with the correct settings"))
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
                setting_string_info(self.rt_unit, "id_file_rt_unit",
                                    settings.id_file_rt_unit, False),
                setting_numerical_info(self.time_window, "time_window",
                                       settings.time_window, False),
                setting_numerical_info(self.ppm_error, "ppm_window",
                                       settings.ppm_window, True),
                setting_string_info(self.use_chromatography_division,
                                    "use_chromatography_division",
                                    settings.use_chromatography_division,
                                    False),
                setting_string_info(self.label_key, "label_key",
                                    settings.label_key, False),
                setting_numerical_info(self.min_AA_length, "min_aa_sequence_length",
                                       settings.min_aa_sequence_length, True),
                setting_numerical_info(self.min_allowed_n_value, "min_allowed_n_values",
                                       settings.min_allowed_n_values, True),
                setting_numerical_info(self.min_allowed_rate, "minimum_allowed_sequence_rate",
                                       settings.minimum_allowed_sequence_rate, False),
                setting_numerical_info(self.max_allowed_rate, "maximum_allowed_sequence_rate",
                                       settings.maximum_allowed_sequence_rate, False),
                setting_string_info(self.graph_file_type,
                                    "graph_output_format",
                                    settings.graph_output_format,
                                    False),
                setting_numerical_info(self.abund_manual_bias, "abundance_manual_bias",
                                       settings.abundance_manual_bias, False),
                setting_string_info(self.asymptope_type, "asymptote",
                                    settings.asymptote, False),
                setting_numerical_info(self.combined_manual_bias, "combined_manual_bias",
                                       settings.combined_manual_bias, False),
                setting_numerical_info(self.fixed_asymptote_value, "fixed_asymptote_value",
                                       settings.fixed_asymptote_value, False),
                setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                       settings.minimum_nonzero_points, True),
                setting_numerical_info(self.proliferation_adjustment, "proliferation_adjustment",
                                       settings.proliferation_adjustment, False),
                setting_string_info(self.use_outlier_removal, "use_outlier_removal",
                                    settings.use_outlier_removal, True),
                setting_numerical_info(self.spacing_manual_bias, "spacing_manual_bias",
                                       settings.spacing_manual_bias, False),
                setting_string_info(self.abundance_type, "abundance_type",
                                    settings.abundance_type, False),
                setting_string_info(self.calculate_n_values, "use_empir_n_value",
                                    settings.use_empir_n_value, True),
                setting_string_info(self.verbose_rate, "verbose_rate",
                                    settings.verbose_rate, True),
                setting_string_info(self.bias_calculation, "bias_calculation",
                                    settings.bias_calculation, False),
                setting_numerical_info(self.max_fn_standard_deviation, "max_fn_standard_deviation",
                                       settings.max_fn_standard_deviation, False),
                setting_string_info(self.fraction_new_calculation, "fraction_new_calculation",
                                    settings.fraction_new_calculation, False),
                setting_numerical_info(self.n_value_cv_limit, "n_value_cv_limit",
                                       settings.n_value_cv_limit, False),
                setting_string_info(self.separate_adducts, "separate_adducts",
                                    settings.separate_adducts, True),
                setting_numerical_info(self.r2_threshold, "r2_threshold",
                                       settings.r2_threshold, False),
                setting_numerical_info(self.zscore_cutoff, "zscore_cutoff",
                                       settings.zscore_cutoff, False),
                setting_string_info(self.graph_n_value_calculations, "graph_n_value_calculations",
                                    settings.graph_n_value_calculations, True),
                setting_string_info(self.save_n_value_data, "save_n_value_data",
                                    settings.save_n_value_data, True)
            ]
            for setting_object in self.all_settings:
                setting_object.set_object_value()
            QtWidgets.QMessageBox.information(self, "Info", "Settings successfully loaded.")
            return

    # load the default peptide settings file
    def load_peptide_default_settings(self):
        settings.load(default_peptide_rate_settings)
        settings.freeze(self.current_setting_file)
        self.all_settings = [
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
            setting_string_info(self.use_chromatography_division,
                                "use_chromatography_division",
                                settings.use_chromatography_division,
                                False),
            setting_string_info(self.label_key, "label_key",
                                settings.label_key, False),
            setting_numerical_info(self.min_AA_length, "min_aa_sequence_length",
                                   settings.min_aa_sequence_length, True),
            setting_numerical_info(self.min_allowed_n_value, "min_allowed_n_values",
                                   settings.min_allowed_n_values, True),
            setting_numerical_info(self.min_allowed_rate, "minimum_allowed_sequence_rate",
                                   settings.minimum_allowed_sequence_rate, False),
            setting_numerical_info(self.max_allowed_rate, "maximum_allowed_sequence_rate",
                                   settings.maximum_allowed_sequence_rate, False),
            setting_string_info(self.graph_file_type,
                                "graph_output_format",
                                settings.graph_output_format,
                                False),
            setting_numerical_info(self.abund_manual_bias, "abundance_manual_bias",
                                   settings.abundance_manual_bias, False),
            setting_string_info(self.asymptope_type, "asymptote",
                                settings.asymptote, False),
            setting_numerical_info(self.combined_manual_bias, "combined_manual_bias",
                                   settings.combined_manual_bias, False),
            setting_numerical_info(self.fixed_asymptote_value, "fixed_asymptote_value",
                                   settings.fixed_asymptote_value, False),
            setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                   settings.minimum_nonzero_points, True),
            setting_numerical_info(self.proliferation_adjustment, "proliferation_adjustment",
                                   settings.proliferation_adjustment, False),
            setting_string_info(self.use_outlier_removal, "use_outlier_removal",
                                settings.use_outlier_removal, True),
            setting_numerical_info(self.spacing_manual_bias, "spacing_manual_bias",
                                   settings.spacing_manual_bias, False),
            setting_string_info(self.abundance_type, "abundance_type",
                                settings.abundance_type, False),
            setting_string_info(self.calculate_n_values, "use_empir_n_value",
                                settings.use_empir_n_value, True),
            setting_string_info(self.verbose_rate, "verbose_rate",
                                settings.verbose_rate, True),
            setting_string_info(self.bias_calculation, "bias_calculation",
                                settings.bias_calculation, False),
            setting_numerical_info(self.max_fn_standard_deviation, "max_fn_standard_deviation",
                                   settings.max_fn_standard_deviation, False),
            setting_string_info(self.fraction_new_calculation, "fraction_new_calculation",
                                settings.fraction_new_calculation, False),
            setting_numerical_info(self.n_value_cv_limit, "n_value_cv_limit",
                                   settings.n_value_cv_limit, False),
            setting_string_info(self.separate_adducts, "separate_adducts",
                                settings.separate_adducts, True),
            setting_numerical_info(self.r2_threshold, "r2_threshold",
                                   settings.r2_threshold, False),
            setting_numerical_info(self.zscore_cutoff, "zscore_cutoff",
                                   settings.zscore_cutoff, False),
            setting_string_info(self.graph_n_value_calculations, "graph_n_value_calculations",
                                settings.graph_n_value_calculations, True),
            setting_string_info(self.save_n_value_data, "save_n_value_data",
                                settings.save_n_value_data, True)
        ]
        for setting_object in self.all_settings:
            setting_object.set_object_value()
        QtWidgets.QMessageBox.information(self, "Info", "Settings successfully loaded.")
        return

    # load the default lipid settings file
    def load_lipid_default_settings(self):
        settings.load(default_lipid_rate_settings)
        settings.freeze(self.current_setting_file)
        self.all_settings = [
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
            setting_string_info(self.use_chromatography_division,
                                "use_chromatography_division",
                                settings.use_chromatography_division,
                                False),
            setting_string_info(self.label_key, "label_key",
                                settings.label_key, False),
            setting_numerical_info(self.min_AA_length, "min_aa_sequence_length",
                                   settings.min_aa_sequence_length, True),
            setting_numerical_info(self.min_allowed_n_value, "min_allowed_n_values",
                                   settings.min_allowed_n_values, True),
            setting_numerical_info(self.min_allowed_rate, "minimum_allowed_sequence_rate",
                                   settings.minimum_allowed_sequence_rate, False),
            setting_numerical_info(self.max_allowed_rate, "maximum_allowed_sequence_rate",
                                   settings.maximum_allowed_sequence_rate, False),
            setting_string_info(self.graph_file_type,
                                "graph_output_format",
                                settings.graph_output_format,
                                False),
            setting_numerical_info(self.abund_manual_bias, "abundance_manual_bias",
                                   settings.abundance_manual_bias, False),
            setting_string_info(self.asymptope_type, "asymptote",
                                settings.asymptote, False),
            setting_numerical_info(self.combined_manual_bias, "combined_manual_bias",
                                   settings.combined_manual_bias, False),
            setting_numerical_info(self.fixed_asymptote_value, "fixed_asymptote_value",
                                   settings.fixed_asymptote_value, False),
            setting_numerical_info(self.minimum_nonzero_points, "minimum_nonzero_points",
                                   settings.minimum_nonzero_points, True),
            setting_numerical_info(self.proliferation_adjustment, "proliferation_adjustment",
                                   settings.proliferation_adjustment, False),
            setting_string_info(self.use_outlier_removal, "use_outlier_removal",
                                settings.use_outlier_removal, True),
            setting_numerical_info(self.spacing_manual_bias, "spacing_manual_bias",
                                   settings.spacing_manual_bias, False),
            setting_string_info(self.abundance_type, "abundance_type",
                                settings.abundance_type, False),
            setting_string_info(self.calculate_n_values, "use_empir_n_value",
                                settings.use_empir_n_value, True),
            setting_string_info(self.verbose_rate, "verbose_rate",
                                settings.verbose_rate, True),
            setting_string_info(self.bias_calculation, "bias_calculation",
                                settings.bias_calculation, False),
            setting_numerical_info(self.max_fn_standard_deviation, "max_fn_standard_deviation",
                                   settings.max_fn_standard_deviation, False),
            setting_string_info(self.fraction_new_calculation, "fraction_new_calculation",
                                settings.fraction_new_calculation, False),
            setting_numerical_info(self.n_value_cv_limit, "n_value_cv_limit",
                                   settings.n_value_cv_limit, False),
            setting_string_info(self.separate_adducts, "separate_adducts",
                                settings.separate_adducts, True),
            setting_numerical_info(self.r2_threshold, "r2_threshold",
                                   settings.r2_threshold, False),
            setting_numerical_info(self.zscore_cutoff, "zscore_cutoff",
                                   settings.zscore_cutoff, False),
            setting_string_info(self.graph_n_value_calculations, "graph_n_value_calculations",
                                settings.graph_n_value_calculations, True),
            setting_string_info(self.save_n_value_data, "save_n_value_data",
                                settings.save_n_value_data, True)
        ]
        for setting_object in self.all_settings:
            setting_object.set_object_value()
        QtWidgets.QMessageBox.information(self, "Info", "Settings successfully loaded.")
        return

    # should overwrite the close of the exit button and the red x in the corner  
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
                allowed_to_save = self.save_settings()
                if allowed_to_save:
                    event.accept()
                else:
                    event.ignore()
                    # hit red x or cancel just don't exit
            else:
                event.ignore()

                # the point of this is to get the values from the settings .yaml

    # that the user is not altering. don't need special classes since
    #  we have to officially declare and all we need to do is get the value
    @staticmethod
    def _get_filters():
        unalterable_settings = {
            "debug_level": settings.debug_level,
            "trim_ids_to_mzml_bounds": settings.trim_ids_to_mzml_bounds,
            "chunk_size": settings.chunk_size,
            "chunking_method_threshold": settings.chunking_method_threshold,
            "peak_ratio_denominator": settings.peak_ratio_denominator,
            "aa_labeling_sites_path": settings.aa_labeling_sites_path,
            "peak_lookback": settings.peak_lookback,
            "peak_lookahead": settings.peak_lookahead,
            "baseline_lookback": settings.baseline_lookback,
            "min_envelopes_to_combine": settings.min_envelopes_to_combine,
            "max_valid_angle": settings.max_valid_angle,
            "max_allowed_enrichment": settings.max_allowed_enrichment,
            "ms_level": settings.ms_level,
            "intensity_filter": settings.intensity_filter,
            "rel_height": settings.rel_height,
            "sampling_rate": settings.sampling_rate,
            "smoothing_width": settings.smoothing_width,
            "smoothing_order": settings.smoothing_order,
            "allowed_peak_variance_min": settings.allowed_peak_variance_min,
            "allowed_neutromer_peak_variance": settings.allowed_neutromer_peak_variance,
            "adduct_weight": settings.adduct_weight,
            "variance_weight": settings.variance_weight,
            "ID_weight": settings.ID_weight,
            "intensity_weight": settings.intensity_weight,
            "how_divided": settings.how_divided,
            "abundance_agreement_filter": settings.abundance_agreement_filter,
            "combined_agreement_filter": settings.combined_agreement_filter,
            "error_of_non_replicated_point": settings.error_of_non_replicated_point,
            "error_of_zero": settings.error_of_zero,
            "lipid_analyte_id_column": settings.lipid_analyte_id_column,
            "lipid_analyte_name_column": settings.lipid_analyte_name_column,
            "peptide_analyte_id_column": settings.peptide_analyte_id_column,
            "peptide_analyte_name_column": settings.peptide_analyte_name_column,
            'remove_filters': settings.remove_filters,
            "s_n_filter": settings.s_n_filter,
            "spacing_agreement_filter": settings.spacing_agreement_filter,
            "unique_sequence_column": settings.unique_sequence_column,
            "y_intercept_of_fit": settings.y_intercept_of_fit
        }
        return unalterable_settings
