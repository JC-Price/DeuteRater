import yaml
import traceback

from pathlib import Path


from utils.exc import InvalidSettingsWarning
from utils.exc import InvalidSettingsError  # noqa: 401 

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

debug_level: int
recognize_available_cores: bool
n_processors: int
id_file_rt_unit: str
trim_ids_to_mzml_bounds: bool
fpe_tolerance: bool
chunk_size: int
chunking_method_threshold: int
max_valid_angle: float
time_window: float
ppm_window: int
heavy_isotope: str
label_key: str
use_abundance: bool
use_neutromer_spacing: bool
maximum_theoretical_pct: float
labeling_step_size: float
minimum_nonzero_points: int
peak_lookback: int
peak_lookahead: int
baseline_lookback: int
#condense_peak_zscore_cutoff: int
peak_ratio_denominator: int
zscore_cutoff: int
mz_proximity_tolerance: float
analyte_id_column: str
analyte_name_column: str
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
use_empir_n_value: bool

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

        global fpe_tolerance
        fpe_tolerance = s['fpe_tolerance']

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

        global heavy_isotope
        heavy_isotope = s['heavy_isotope']

        global label_key
        label_key = s["label_key"]

        global use_abundance
        use_abundance = s['use_abundance']

        global use_neutromer_spacing
        use_neutromer_spacing = s['use_neutromer_spacing']

        global maximum_theoretical_pct
        maximum_theoretical_pct = s['maximum_theoretical_pct']

        global labeling_step_size
        labeling_step_size = s['labeling_step_size']

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

        global mz_proximity_tolerance
        mz_proximity_tolerance = s['mz_proximity_tolerance']
        
        global analyte_id_column
        analyte_id_column = s['analyte_id_column']
        
        global analyte_name_column
        analyte_name_column = s['analyte_name_column']
        
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

        global use_empir_n_value
        use_empir_n_value = s["use_empir_n_value"]
        

    except Exception as e:
        print(e)
        traceback.print_tb(e.__traceback__)


def freeze(path=None, settings_dict = None):
    if not settings_dict:
        settings_dict = {
            'debug_level': debug_level,
            'recognize_available_cores': recognize_available_cores,
            'n_processors': n_processors,
            'id_file_rt_unit': id_file_rt_unit,
            'trim_ids_to_mzml_bounds': trim_ids_to_mzml_bounds,
            'fpe_tolerance': fpe_tolerance,
            'chunk_size': chunk_size,
            'chunking_method_threshold': chunking_method_threshold,
            'max_valid_angle': max_valid_angle,
            'time_window': time_window,
            'ppm_window': ppm_window,
            'heavy_isotope': heavy_isotope,
            "label_key": label_key,
            'use_abundance': use_abundance,
            'use_neutromer_spacing': use_neutromer_spacing,
            'maximum_theoretical_pct': maximum_theoretical_pct,
            'labeling_step_size': labeling_step_size,
            'minimum_nonzero_points': minimum_nonzero_points,
            'peak_lookback': peak_lookback,
            'peak_lookahead': peak_lookahead,
            'baseline_lookback': baseline_lookback,
            'min_envelopes_to_combine': min_envelopes_to_combine,
            'peak_ratio_denominator': peak_ratio_denominator,
            'zscore_cutoff': zscore_cutoff,
            'mz_proximity_tolerance': mz_proximity_tolerance,
            "analyte_id_column": analyte_id_column,
            "analyte_name_column": analyte_name_column,
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
            "min_allowed_n_values": min_allowed_n_values
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

