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


import yaml
import traceback
import os

from pathlib import Path

try:
    from utils.exc import InvalidSettingsWarning
    from utils.exc import InvalidSettingsError  # noqa: 401
except:
    from DeuteRater.utils.exc import InvalidSettingsError, InvalidSettingsWarning

# TODO: How would I dynamically load a different settings file?
# TODO: ^^^This really should be able to be passed in
# TODO: should different steps have different settings files?
# TODO: add reasonable constraints on settings
# TODO: verbose exception output
# TODO: add documentation on what this file is for
# TODO: Shorten variable names where possible
# TODO: add error checking where applicable
# TODO: set type annotations and default values
# TODO: determine good defaults
# TODO: type annotations for path variable
# TODO: Discuss this style vs settings class style
# TODO: Figure out how to reduce redundancy. Like a better singleton

# NOTE: we can use LibYAML c bindings if we need more speed

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
resource_location = os.path.join(location, "resources")

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
study_type: str
aa_label_path: str
use_abundance: str
use_neutromer_spacing: bool
minimum_nonzero_points: int
peak_lookback: int
peak_lookahead: int
baseline_lookback: int
peak_ratio_denominator: int
zscore_cutoff: int
rt_proximity_tolerance: float
mz_proximity_tolerance: float
peptide_analyte_id_column: str
peptide_analyte_name_column: str
unique_sequence_column: str
roll_up_rate_calc: bool
asymptote: str
proliferation_adjustment: float
fixed_asymptote_value: float
error_of_zero: float
abundance_agreement_filter: float
spacing_agreement_filter: float
combined_agreement_filter: float
bias_calculation: str
abundance_manual_bias: float
spacing_manual_bias: float
combined_manual_bias: float
y_intercept_of_fit: float
error_of_non_replicated_point: float
enrichement_of_zero: float
min_allowed_abund_max_delta: float
min_aa_sequence_length: int
min_allowed_n_values: int
minimum_abund_change: float
intensity_filter: int
verbose_rate: bool
rate_output_format: str


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
        
        global study_type
        study_type = s['study_type']
        
        global aa_label_path
        aa_label_path = os.path.join(resource_location,
                                     s['aa_labeling_sites_path'])

        global use_abundance
        use_abundance = s['use_abundance']
        
        global use_neutromer_spacing
        use_neutromer_spacing = s['use_neutromer_spacing']
        
        global minimum_nonzero_points
        minimum_nonzero_points = s['minimum_nonzero_points']
        
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
        
        
        global rt_proximity_tolerance
        rt_proximity_tolerance = s['rt_proximity_tolerance']
        
        global mz_proximity_tolerance
        mz_proximity_tolerance = s['mz_proximity_tolerance']
        
        global peptide_analyte_id_column
        peptide_analyte_id_column = s['peptide_analyte_id_column']
        
        global peptide_analyte_name_column
        peptide_analyte_name_column = s['peptide_analyte_name_column']
        
        global unique_sequence_column
        unique_sequence_column = s["unique_sequence_column"]
        
        global roll_up_rate_calc
        roll_up_rate_calc = s['roll_up_rate_calc']
        
        global asymptote
        asymptote = s['asymptote']
        
        global proliferation_adjustment
        proliferation_adjustment = s['proliferation_adjustment']
        
        global fixed_asymptote_value
        fixed_asymptote_value = s["fixed_asymptote_value"]
        
        global error_of_zero
        error_of_zero = s["error_of_zero"]
        
        global abundance_agreement_filter
        abundance_agreement_filter = s["abundance_agreement_filter"]
        
        global spacing_agreement_filter
        spacing_agreement_filter = s["spacing_agreement_filter"]
        
        global combined_agreement_filter
        combined_agreement_filter = s["combined_agreement_filter"]
        
        global bias_calculation
        bias_calculation = s["bias_calculation"]
        
        global abundance_manual_bias
        abundance_manual_bias = s["abundance_manual_bias"]
        
        global spacing_manual_bias
        spacing_manual_bias = s["spacing_manual_bias"]
        
        global combined_manual_bias
        combined_manual_bias = s["combined_manual_bias"]
        
        global y_intercept_of_fit
        y_intercept_of_fit = s["y_intercept_of_fit"]
        
        global error_of_non_replicated_point
        error_of_non_replicated_point = s["error_of_non_replicated_point"]
        
        global enrichement_of_zero
        enrichement_of_zero = s["enrichement_of_zero"]
        
        global min_allowed_abund_max_delta
        min_allowed_abund_max_delta = s["min_allowed_abund_max_delta"]
        
        global min_aa_sequence_length
        min_aa_sequence_length = s["min_aa_sequence_length"]
        
        global min_allowed_n_values
        min_allowed_n_values = s["min_allowed_n_values"]
        
        global minimum_abund_change
        minimum_abund_change = s["minimum_abund_change"]
        
        global intensity_filter
        intensity_filter = s["intensity_filter"]

        global verbose_rate
        verbose_rate = s["verbose_rate"]
        
        global ms_level
        ms_level = s["ms_level"]
        
        global rate_output_format
        rate_output_format = s["rate_output_format"]
        
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


def freeze(path=None, settings_dict=None):
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
            'study_type': study_type,
            'aa_labeling_sites_path': aa_label_path,
            'use_abundance': use_abundance,
            'use_neutromer_spacing': use_neutromer_spacing,
            'minimum_nonzero_points': minimum_nonzero_points,
            'peak_lookback': peak_lookback,
            'peak_lookahead': peak_lookahead,
            'baseline_lookback': baseline_lookback,
            'min_envelopes_to_combine': min_envelopes_to_combine,
            'peak_ratio_denominator': peak_ratio_denominator,
            'zscore_cutoff': zscore_cutoff,
            "rt_proximity_tolerance": rt_proximity_tolerance,
            'mz_proximity_tolerance': mz_proximity_tolerance,
            "peptide_analyte_id_column": peptide_analyte_id_column,
            "peptide_analyte_name_column": peptide_analyte_name_column,
            "unique_sequence_column": unique_sequence_column,
            "roll_up_rate_calc": roll_up_rate_calc,
            "asymptote": asymptote,
            "proliferation_adjustment": proliferation_adjustment,
            "fixed_asymptote_value": fixed_asymptote_value,
            "error_of_zero": error_of_zero,
            "abundance_agreement_filter": abundance_agreement_filter,
            "spacing_agreement_filter": spacing_agreement_filter,
            "combined_agreement_filter": combined_agreement_filter,
            "bias_calculation": bias_calculation,
            "abundance_manual_bias": abundance_manual_bias,
            "spacing_manual_bias": spacing_manual_bias,
            "combined_manual_bias": combined_manual_bias,
            "y_intercept_of_fit": y_intercept_of_fit,
            "error_of_non_replicated_point": error_of_non_replicated_point,
            "enrichement_of_zero": enrichement_of_zero,
            "min_allowed_abund_max_delta": min_allowed_abund_max_delta,
            "min_aa_sequence_length": min_aa_sequence_length,
            "min_allowed_n_values": min_allowed_n_values,
            "minimum_abund_change": minimum_abund_change,
            "verbose_rate": verbose_rate,
            "intensity_filter": intensity_filter,
            "ms_level": ms_level,
            "rate_output_format": rate_output_format
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
