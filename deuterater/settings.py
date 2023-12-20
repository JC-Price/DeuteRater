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
starting_enrichment_table_timepoints: int
error_estimation: str
min_non_zero_timepoints_rate: int
min_allowed_timepoints_enrichment: int
minimum_allowed_sequence_rate: float
maximum_allowed_sequence_rate: float
minimum_sequences_to_combine_for_protein_rate: int
lowest_allowed_norm_isotope: float
highest_allowed_norm_isotope: float
m0_decreasing_allowed_noise: float
median_absolute_residuals_cutoff_single_point: float
median_absolute_residuals_cutoff_two_points: float
median_absolute_residuals_cutoff_general: float
desired_points_for_optimization_graph: int
intensity_filter: int
rel_height: float
protein_combination_method: str
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
max_enrichment_allowed: float
verbose_output: bool

# TODO: add quick explanation of how this works, inc. 'global' doc link
def load(settings_path):
    # NOTE: Look at the python documentation for the 'global' statement if youf
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
        
        global starting_enrichment_table_timepoints
        starting_enrichment_table_timepoints = s["starting_enrichment_table_timepoints"]
        
        global error_estimation
        error_estimation = s["error_estimation"]
        
        global min_non_zero_timepoints_rate
        min_non_zero_timepoints_rate = s["min_non_zero_timepoints_rate"]
        
        global min_allowed_timepoints_enrichment
        min_allowed_timepoints_enrichment = s["min_allowed_timepoints_enrichment"]
        
        global minimum_allowed_sequence_rate
        minimum_allowed_sequence_rate = s["minimum_allowed_sequence_rate"]
        
        global maximum_allowed_sequence_rate
        maximum_allowed_sequence_rate = s["maximum_allowed_sequence_rate"]
        
        global minimum_sequences_to_combine_for_protein_rate
        minimum_sequences_to_combine_for_protein_rate = s["minimum_sequences_to_combine_for_protein_rate"]
        
        global lowest_allowed_norm_isotope
        lowest_allowed_norm_isotope = s["lowest_allowed_norm_isotope"]
        
        global highest_allowed_norm_isotope
        highest_allowed_norm_isotope = s["highest_allowed_norm_isotope"]
        
        global m0_decreasing_allowed_noise
        m0_decreasing_allowed_noise = s["m0_decreasing_allowed_noise"]
        
        global median_absolute_residuals_cutoff_single_point
        median_absolute_residuals_cutoff_single_point = s["median_absolute_residuals_cutoff_single_point"]
        
        global median_absolute_residuals_cutoff_two_points
        median_absolute_residuals_cutoff_two_points = s["median_absolute_residuals_cutoff_two_points"]
        
        global median_absolute_residuals_cutoff_general
        median_absolute_residuals_cutoff_general = s["median_absolute_residuals_cutoff_general"]
        
        global desired_points_for_optimization_graph
        desired_points_for_optimization_graph = s["desired_points_for_optimization_graph"]
        
        global intensity_filter
        intensity_filter = s["intensity_filter"]
        
        global rel_height
        rel_height = s["rel_height"]
        
        global sampling_rate
        sampling_rate = s["sampling_rate"]
        
        global protein_combination_method
        protein_combination_method = s["protein_combination_method"]
        
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
        
        global max_enrichment_allowed
        max_enrichment_allowed = s["max_enrichment_allowed"]
        
        global verbose_output
        verbose_output = s["verbose_output"]

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
            "max_enrichment_allowed": max_enrichment_allowed,
            "min_aa_sequence_length": min_aa_sequence_length,
            "mz_proximity_tolerance":mz_proximity_tolerance,
            "rt_proximity_tolerance":rt_proximity_tolerance,
            "min_allowed_n_values": min_allowed_n_values,
            "starting_enrichment_table_timepoints": starting_enrichment_table_timepoints,
            "error_estimation": error_estimation,
            "min_non_zero_timepoints_rate": min_non_zero_timepoints_rate,
            "min_allowed_timepoints_enrichment": min_allowed_timepoints_enrichment,
            "minimum_allowed_sequence_rate": minimum_allowed_sequence_rate,
            "maximum_allowed_sequence_rate": maximum_allowed_sequence_rate,
            "minimum_sequences_to_combine_for_protein_rate": minimum_sequences_to_combine_for_protein_rate,
            "lowest_allowed_norm_isotope": lowest_allowed_norm_isotope,
            "highest_allowed_norm_isotope": highest_allowed_norm_isotope,
            "m0_decreasing_allowed_noise": m0_decreasing_allowed_noise,
            "median_absolute_residuals_cutoff_single_point": median_absolute_residuals_cutoff_single_point,
            "median_absolute_residuals_cutoff_two_points": median_absolute_residuals_cutoff_two_points,
            "median_absolute_residuals_cutoff_general": median_absolute_residuals_cutoff_general,
            "desired_points_for_optimization_graph": desired_points_for_optimization_graph,
            "intensity_filter": intensity_filter,
            "rel_height": rel_height,
            "sampling_rate": sampling_rate,
            "protein_combination_method": protein_combination_method,
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
            "verbose_output": verbose_output
            
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

