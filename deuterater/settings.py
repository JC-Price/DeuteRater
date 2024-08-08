# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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


import yaml
import traceback
import os

from pathlib import Path


from utils.exc import InvalidSettingsWarning
from utils.exc import InvalidSettingsError  # noqa: 401 

"""
this is to load the settings as global variables.  There are three portions
the first defines the variable type
the second is the load function which loads the variables
the third is the freeze function which allows saving the variables

"""

location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
resource_location = os.path.join(location, "resources")

# TODO: go through settings and make sure they can be used for both proteins and lipids. - Ben D
debug_level: int
recognize_available_cores: bool
n_processors: int
id_file_rt_unit: str
trim_ids_to_mzml_bounds: bool
chunk_size: int
chunking_method_threshold: int
max_valid_angle: float
time_window: float
ppm_window: int
label_key: str
aa_labeling_sites_path: str
peak_lookback: int
peak_lookahead: int
baseline_lookback: int
min_envelopes_to_combine: int
peak_ratio_denominator: int
zscore_cutoff: int
mz_proximity_tolerance: int
rt_proximity_tolerance: float
min_aa_sequence_length: int
min_allowed_n_values: int
minimum_allowed_sequence_rate: float
maximum_allowed_sequence_rate: float
intensity_filter: int
rel_height: float
sampling_rate: int
smoothing_width: int
smoothing_order: int
allowed_peak_variance_min: float
allowed_neutromer_peak_variance: float
adduct_weight: float
variance_weight: float
ID_weight: float
intensity_weight: float
how_divided: str
use_chromatography_division: str
graph_output_format: str
ms_level: int
abundance_agreement_filter: float
abundance_manual_bias: float
asymptote: str
combined_agreement_filter: float
combined_manual_bias: float
error_of_non_replicated_point: float
error_of_zero: float
fixed_asymptote_value: float
lipid_analyte_id_column: str
lipid_analyte_name_column: str
minimum_nonzero_points: int
peptide_analyte_id_column: str
peptide_analyte_name_column: str
proliferation_adjustment: float
remove_filters: bool
roll_up_rate_calc: bool
s_n_filter: float
separate_adducts: bool
spacing_agreement_filter: float
spacing_manual_bias: float
unique_sequence_column: str
use_empir_n_value: bool
abundance_type: str
bias_calculation: str
verbose_rate: bool
y_intercept_of_fit: float
max_fn_standard_deviation: float
fraction_new_calculation: str


# TODO: add quick explanation of how this works, inc. 'global' doc link
def load(settings_path):
    # NOTE: Look at the python documentation for the 'global' statement if you
    #       Want to understand how this module works
    try:
        settings_path = Path(settings_path)
        with settings_path.open('r') as f:
            s = yaml.load(f, Loader=yaml.FullLoader)
        global debug_level
        debug_level = s['debug_level']
        if debug_level not in [0, 1, 2]:
            raise InvalidSettingsWarning(
                'Invalid debug level value given'
            )
            print('Running with debug level 0')
            debug_level = 0

        global recognize_available_cores
        recognize_available_cores = s['recognize_available_cores']

        global n_processors
        n_processors = s['n_processors']

        global id_file_rt_unit
        id_file_rt_unit = s['id_file_rt_unit']

        global trim_ids_to_mzml_bounds
        trim_ids_to_mzml_bounds = s['trim_ids_to_mzml_bounds']

        global chunk_size
        chunk_size = s['chunk_size']

        global chunking_method_threshold
        chunking_method_threshold = s['chunking_method_threshold']
        
        global max_valid_angle
        max_valid_angle = s['max_valid_angle']

        global time_window
        time_window = s['time_window']

        global ppm_window
        ppm_window = s['ppm_window']

        global label_key
        label_key = s["label_key"]
        
        global aa_labeling_sites_path
        aa_labeling_sites_path = os.path.join(resource_location, s["aa_labeling_sites_path"])

        global peak_lookback
        peak_lookback = s['peak_lookback']

        global peak_lookahead
        peak_lookahead = s['peak_lookahead']

        global baseline_lookback
        baseline_lookback = s['baseline_lookback']

        global min_envelopes_to_combine
        min_envelopes_to_combine = s['min_envelopes_to_combine']

        global peak_ratio_denominator
        peak_ratio_denominator = s['peak_ratio_denominator']

        global zscore_cutoff
        zscore_cutoff = s['zscore_cutoff']
        
        global mz_proximity_tolerance
        mz_proximity_tolerance = s["mz_proximity_tolerance"]
        
        global rt_proximity_tolerance
        rt_proximity_tolerance = s["rt_proximity_tolerance"]
        
        global min_aa_sequence_length
        min_aa_sequence_length = s["min_aa_sequence_length"]
        
        global min_allowed_n_values
        min_allowed_n_values = s["min_allowed_n_values"]
        
        global minimum_allowed_sequence_rate
        minimum_allowed_sequence_rate = s["minimum_allowed_sequence_rate"]
        
        global maximum_allowed_sequence_rate
        maximum_allowed_sequence_rate = s["maximum_allowed_sequence_rate"]

        global intensity_filter
        intensity_filter = s["intensity_filter"]
        
        global rel_height
        rel_height = s["rel_height"]
        
        global sampling_rate
        sampling_rate = s["sampling_rate"]
        
        global smoothing_width
        smoothing_width = s["smoothing_width"]
        
        global smoothing_order
        smoothing_order = s["smoothing_order"]
        
        global allowed_peak_variance_min
        allowed_peak_variance_min = s["allowed_peak_variance_min"]
        
        global adduct_weight
        adduct_weight = s["adduct_weight"]
        
        global variance_weight
        variance_weight = s["variance_weight"]
        
        global ID_weight
        ID_weight = s["ID_weight"]
        
        global intensity_weight
        intensity_weight = s["intensity_weight"]
        
        global how_divided
        how_divided = s["how_divided"]
        
        global allowed_neutromer_peak_variance
        allowed_neutromer_peak_variance = s["allowed_neutromer_peak_variance"]
        
        global ms_level
        ms_level = s["ms_level"]
        
        global use_chromatography_division
        use_chromatography_division = s["use_chromatography_division"]
        
        global graph_output_format
        graph_output_format = s["graph_output_format"]

        global abundance_agreement_filter
        abundance_agreement_filter = s["abundance_agreement_filter"]

        global abundance_manual_bias
        abundance_manual_bias = s["abundance_manual_bias"]

        global asymptote
        asymptote = s['asymptote']

        global combined_agreement_filter
        combined_agreement_filter = s["combined_agreement_filter"]

        global combined_manual_bias
        combined_manual_bias = s["combined_manual_bias"]

        global error_of_non_replicated_point
        error_of_non_replicated_point = s["error_of_non_replicated_point"]

        global error_of_zero
        error_of_zero = s["error_of_zero"]

        global fixed_asymptote_value
        fixed_asymptote_value = s["fixed_asymptote_value"]

        global lipid_analyte_id_column
        lipid_analyte_id_column = s['lipid_analyte_id_column']

        global lipid_analyte_name_column
        lipid_analyte_name_column = s['lipid_analyte_name_column']

        global minimum_nonzero_points
        minimum_nonzero_points = s['minimum_nonzero_points']

        global peptide_analyte_id_column
        peptide_analyte_id_column = s['peptide_analyte_id_column']

        global peptide_analyte_name_column
        peptide_analyte_name_column = s['peptide_analyte_name_column']

        global proliferation_adjustment
        proliferation_adjustment = s['proliferation_adjustment']

        global remove_filters
        remove_filters = s["remove_filters"]

        global roll_up_rate_calc
        roll_up_rate_calc = s['roll_up_rate_calc']

        global s_n_filter
        s_n_filter = s["s_n_filter"]

        global separate_adducts
        separate_adducts = s["separate_adducts"]

        global spacing_agreement_filter
        spacing_agreement_filter = s["spacing_agreement_filter"]

        global spacing_manual_bias
        spacing_manual_bias = s["spacing_manual_bias"]

        global unique_sequence_column
        unique_sequence_column = s["unique_sequence_column"]

        global use_empir_n_value
        use_empir_n_value = s["use_empir_n_value"]

        global abundance_type
        abundance_type = s['abundance_type']

        global bias_calculation
        bias_calculation = s['bias_calculation']

        global verbose_rate
        verbose_rate = s["verbose_rate"]

        global y_intercept_of_fit
        y_intercept_of_fit = s["y_intercept_of_fit"]

        global max_fn_standard_deviation
        max_fn_standard_deviation = s["max_fn_standard_deviation"]

        global fraction_new_calculation
        fraction_new_calculation = s["fraction_new_calculation"]

    except Exception as e:
        print(e)
        traceback.print_tb(e.__traceback__)


def compare(settings_path, compare_path):
    try:
        settings_path = Path(settings_path)
        with settings_path.open('r') as f:
            setting = yaml.load(f, Loader=yaml.FullLoader)
        compare_path = Path(compare_path)
        with compare_path.open('r') as f:
            compare = yaml.load(f, Loader=yaml.FullLoader)
        if setting.keys() != compare.keys():
            return "Different Keys"
        for key in setting.keys():
            if setting[key] != compare[key]:
                return "Mismatched Keys"
        return "MATCH"
    except:
        return "Error"


def freeze(path=None, settings_dict = None):
    if not settings_dict:
        settings_dict = {
            'debug_level': debug_level,
            'recognize_available_cores': recognize_available_cores,
            'n_processors': n_processors,
            'id_file_rt_unit': id_file_rt_unit,
            'trim_ids_to_mzml_bounds': trim_ids_to_mzml_bounds,
            'chunk_size': chunk_size,
            'chunking_method_threshold': chunking_method_threshold,
            'max_valid_angle': max_valid_angle,
            'time_window': time_window,
            'ppm_window': ppm_window,
            "label_key": label_key,
            "aa_labeling_sites_path": aa_labeling_sites_path,
            'peak_lookback': peak_lookback,
            'peak_lookahead': peak_lookahead,
            'baseline_lookback': baseline_lookback,
            'min_envelopes_to_combine': min_envelopes_to_combine,
            'peak_ratio_denominator': peak_ratio_denominator,
            'zscore_cutoff': zscore_cutoff,
            "min_aa_sequence_length": min_aa_sequence_length,
            "mz_proximity_tolerance":mz_proximity_tolerance,
            "rt_proximity_tolerance":rt_proximity_tolerance,
            "min_allowed_n_values": min_allowed_n_values,
            "minimum_allowed_sequence_rate": minimum_allowed_sequence_rate,
            "maximum_allowed_sequence_rate": maximum_allowed_sequence_rate,
            "intensity_filter": intensity_filter,
            "rel_height": rel_height,
            "sampling_rate": sampling_rate,
            "smoothing_width": smoothing_width,
            "smoothing_order": smoothing_order,
            "allowed_peak_variance_min": allowed_peak_variance_min,
            "adduct_weight": adduct_weight,
            "variance_weight": variance_weight,
            "ID_weight": ID_weight,
            "intensity_weight": intensity_weight,
            "how_divided": how_divided,
            "allowed_neutromer_peak_variance": allowed_neutromer_peak_variance,
            "ms_level": ms_level,
            "use_chromatography_division": use_chromatography_division,
            "graph_output_format": graph_output_format,
            "abundance_agreement_filter": abundance_agreement_filter,
            "abundance_manual_bias": abundance_manual_bias,
            "asymptote": asymptote,
            "combined_agreement_filter": combined_agreement_filter,
            "combined_manual_bias": combined_manual_bias,
            "error_of_non_replicated_point": error_of_non_replicated_point,
            "error_of_zero": error_of_zero,
            "fixed_asymptote_value": fixed_asymptote_value,
            "lipid_analyte_id_column": lipid_analyte_id_column,
            "lipid_analyte_name_column": lipid_analyte_name_column,
            "minimum_nonzero_points": minimum_nonzero_points,
            "peptide_analyte_id_column": peptide_analyte_id_column,
            "peptide_analyte_name_column": peptide_analyte_name_column,
            "proliferation_adjustment": proliferation_adjustment,
            "remove_filters": remove_filters,
            "roll_up_rate_calc": roll_up_rate_calc,
            "s_n_filter": s_n_filter,
            "separate_adducts": separate_adducts,
            "spacing_agreement_filter": spacing_agreement_filter,
            "spacing_manual_bias": spacing_manual_bias,
            "unique_sequence_column": unique_sequence_column,
            "use_empir_n_value": use_empir_n_value,
            "abundance_type": abundance_type,
            "bias_calculation": bias_calculation,
            "verbose_rate": verbose_rate,
            "y_intercept_of_fit": y_intercept_of_fit,
            "max_fn_standard_deviation": max_fn_standard_deviation,
            "fraction_new_calculation": fraction_new_calculation
        }
    if path:
        with open(path, 'w') as frozen_settings_file:
            yaml.dump(
                data=settings_dict,
                stream=frozen_settings_file,
                canonical=False
            )
    else:
        print(yaml.dump(data=settings_dict, canonical=False))

