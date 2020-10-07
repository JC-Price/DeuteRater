import pymzml
import warnings
import pandas as pd
from numpy import median

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.mzml as dml
import utils.math as dmt
from utils.exc import EmptyIdChunkWarning
from constants import NEUTRON

from obs.peak import Peak
from obs.envelope import Envelope
from obs.id import ID

# TODO: Any temporary values need to be in settings
# TODO: abstract some stepts into functions
# TODO: Research how to properly handle exceptions when multiprocessing
# TODO: stress what 'spectrum_index' and 'native_id' mean
# TODO: should we delete invalid envelopes?

# NOTE: this code assumes that an mzml is ordered by retention time
# TODO: should we add logic to check if the mzml is time ordered?


def extract(settings_path, mzml_path, index_to_ID, chunk):
    '''Extract data from the mzml according to the identification information
    '''
    # A rough outline of how the logic flows.
    # EXTRACT:
    # for each scan in the reader
    #   in ms level 1
    #   for each id in the window
    #       isotope extraction specific logic

    # Turn the warnings off so that it doesn't mess up the tqdm progress bar
    warnings.filterwarnings("ignore")

    # Load the settings into the settings module. This needs to be done in each
    #   process due to how python handles multiprocessing
    settings.load(settings_path)

    # Open a file pointer to the mzml file
    mzml_fp = pymzml.run.Reader(
        path_or_file=mzml_path,
        build_index_from_scratch=True
    )
    mzml_bounds = dml.get_bounds(mzml_fp, index_to_ID)

    # Check for an empty id file chunk
    if len(chunk) <= 0:
        warnings.warn(EmptyIdChunkWarning(
            'There are no identifications in this chunk'
        ))

    # Set the high and low retention time bounds, based on the chunk of the
    #   identification file
    lo_rt_bound = chunk.at[0, 'rt'] - settings.time_window
    if lo_rt_bound < 0:
        lo_rt_bound = 0
    hi_rt_bound = chunk.at[len(chunk)-1, 'rt'] + settings.time_window

    # Search for the scans at the high and low retention time bounds
    lo_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, lo_rt_bound)
    hi_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, hi_rt_bound)
##    if mzml_fp[index_to_ID[hi_spec_idx]] > hi_rt_bound:
##        hi_spec_idx = hi_spec_idx - 1

    # Logic block for handling out-of-bounds indices
    if lo_spec_idx != -1 and hi_spec_idx != -1:
        # Do nothing if both indices are in bounds
        pass
    elif lo_spec_idx == -1 and hi_spec_idx != -1:
        # If the just the higher index is found, assign the lowest index in
        #   the mzml to 'lo_spec_idx'
        lo_spec_idx = mzml_bounds['idx_min']
    elif lo_spec_idx != -1 and hi_spec_idx == -1:
        # If the just the lower index is found, assign the highest index in
        #   the mzml to 'hi_spec_idx'
        hi_spec_idx = mzml_bounds['idx_max']
    elif lo_rt_bound < mzml_bounds['rt_min'] < \
            mzml_bounds['rt_max'] < hi_rt_bound:
        # If neither index is found but the time span covered by the chunk of
        #   the ID file encompasses that of the mzml, assign 'lo_spec_idx' and
        #   'hi_spec_idx' the minimum and maximum index values given by the
        #   mzml file
        lo_spec_idx = mzml_bounds['idx_min']
        hi_spec_idx = mzml_bounds['idx_max']
    else:
        # Otherwise, there is no intersection between the ID file and the mzml
        #   in terms of retention time and no analysis can be made
        return -1

    ids = []  # initialize the list of identifications

    # TODO: redefine this column as ionmass?
    chunk['mass'] = chunk['mz'] * chunk['z']

    # Instantiate all of the identifications in the chunk
    for row in chunk.itertuples(index=True):
        ids.append(
            ID(
                rt=row.rt,
                mz=row.mz,
                mass=row.mass,
                z=row.z,
                n_isos=row.n_isos,
                #cf=row.cf
            )
        )

    # Iterate through all of the relevent spectrums in the mzml
    for spectrum_index in dmt.inclusive_range(lo_spec_idx, hi_spec_idx):
        # apply the index_to_ID map in order to access the correct spectrum
        native_id = index_to_ID[spectrum_index]
        try:
            # try to access this spectrum
            spectrum = mzml_fp[native_id]
            spec_rt = spectrum.scan_time_in_minutes()
            spec_mzs = spectrum.mz
            spec_abs = spectrum.i
        except Exception:
            # TODO: use a more specific Exception
            # catch the exception and move on if the spectru is not found
            continue

        # only deal with ms_level 1 for now
        if spectrum.ms_level != 1:
            continue

        # determine id indices of peaks searches
        # adding and subtracting the floating point error tolerance allows us
        # to include the extremes of the range
        local_window_min = \
            spec_rt - (settings.time_window)  # + settings.fpe_tolerance)
        local_window_max = \
            spec_rt + (settings.time_window)  # + settings.fpe_tolerance)
        try:
            lo_slice_index = \
                min(chunk[chunk['rt'] > local_window_min].axes[0].tolist())
            hi_slice_index = \
                max(chunk[chunk['rt'] < local_window_max].axes[0].tolist())
        except:
            continue

        # iterate through relevant ids
        for id in ids[dmt.inclusive_slice(lo_slice_index, hi_slice_index)]:
            charge = id.z
            # instatntiate an envelope
            envelope = Envelope(
                peaks=[],
                rt=spec_rt,
                n_lookback=settings.peak_lookback,
                n_lookahead=settings.peak_lookahead
            )

            lo_baseline_bound = None
            hi_baseline_bound = None

            peak_range_start = 0 - settings.peak_lookback
            peak_range_end = id.n_isos + settings.peak_lookahead

            # Iterate through all of the peaks we want to look for
            for peak_num in range(peak_range_start, peak_range_end):
                # define the mz to search for in the spectrum
                search_mz = id.mz + (peak_num * NEUTRON / charge)
                # define the ppm error tolerance
                reach = settings.ppm_window / 1_000_000.0 * search_mz
                # find the index of the nearest data point in that spectrum's
                #   mz array
                index = dmt.find_nearest_index(spec_mzs, search_mz)

                if peak_num == 0:
                    # set the bounds for defining the baseline
                    lo_baseline_bound = dmt.find_nearest_index(
                        spec_mzs,
                        id.mz - settings.baseline_lookback
                    )
                    hi_baseline_bound = index

                # TODO: Do I need to speed this up by removing typecheck?
                # TODO: Expand this to only one paren/bracket per line?
                if abs(spec_mzs[index] - search_mz) < reach:
                    # If the value at that index is within the reach
                    envelope.append_peak(Peak(
                        mz=spec_mzs[index],
                        abundance=spec_abs[index],
                        i=peak_num
                    ))
                else:
                    if 0 <= peak_num < id.n_isos:
                        # set the envelopes validity flag to false if no peak
                        #   is found, then move on to the next identification
                        envelope.is_valid = False
                        break
                    else:
                        # Unless it is one of the extra peaks
                        envelope.append_peak(Peak(
                            mz=search_mz,
                            # TODO: it might be better to set this to NA
                            abundance=0,
                            i=peak_num
                        ))

            # TODO Do i need to speed this up by removing typecheck?
            if envelope.is_valid:
                # If all of the peaks have been found, add it to the
                #   identification (after determining the baseline)
                # NOTE: baseline is defined as the median abundance of the 100
                #   mz units preceding the m0 peak
                envelope.baseline = median(spec_mzs[
                    dmt.inclusive_slice(lo_baseline_bound, hi_baseline_bound)
                ])
                id.append_envelope(envelope)
    mzml_fp.close()

    for id in ids:
        id.aggregate_envelopes()

    # TODO: Better variable naming here. obs? I can do better
    # TODO: is there better way to initialize this?
    # TODO: add lookback columns?

    # Initialize the dataframe to send back to the main process
    peak_obs = pd.DataFrame(
        index=chunk.index.values,
        columns=['mzs', 'abundances',
                 'lookback_mzs', 'lookback_abundances',
                 'lookahead_mzs', 'lookahead_abundances',
                 'rt_min', 'rt_max', 'baseline_signal', 'mads',
                 # 'envelopes_before_angle_filter',
                 # 'envelopes_after_angle_filter',
                 # 'envelopes_outside_angle_filter',
                 # 'deviation_before_angle_filter',
                 # 'deviation_after_angle_filter',
                 # 'deviation_outside_angle_filter',
                 # 'max_m0_abundance',
                 'mzml_path']
    )

    # Populate the
    for row in peak_obs.itertuples():
        i = row.Index
        # TODO: rename condensed envelope to found envelope?
        if ids[i].condensed_envelope:
            # this will not run if condensed_envelope is still 'None'
            mzs, abundances = ids[i].condensed_envelope.to_obs()
            lb_mzs, lb_abundances = ids[i].condensed_envelope.lb_obs()
            la_mzs, la_abundances = ids[i].condensed_envelope.la_obs()
            peak_obs.at[i, 'mzs'] = mzs
            peak_obs.at[i, 'abundances'] = abundances
            peak_obs.at[i, 'lookback_mzs'] = lb_mzs
            peak_obs.at[i, 'lookback_abundances'] = lb_abundances
            peak_obs.at[i, 'lookahead_mzs'] = la_mzs
            peak_obs.at[i, 'lookahead_abundances'] = la_abundances
            peak_obs.at[i, 'rt_min'] = ids[i].rt_min
            peak_obs.at[i, 'rt_max'] = ids[i].rt_max
            peak_obs.at[i, 'baseline_signal'] = \
                ids[i].condensed_envelope.baseline
            peak_obs.at[i, 'mads'] = tuple(ids[i].mads)
            # peak_obs.at[i, 'envelopes_before_angle_filter'] = \
            #     ids[i].envelopes_before_angle_filter
            # peak_obs.at[i, 'envelopes_after_angle_filter'] = \
            #     ids[i].envelopes_after_angle_filter
            # Vectorized math removed from for loop, see below
            # peak_obs.at[i, 'deviation_before_angle_filter'] = \
            #     ids[i].deviation_before_angle_filter
            # peak_obs.at[i, 'deviation_after_angle_filter'] = \
            #     ids[i].deviation_after_angle_filter
            # peak_obs.at[i, 'deviation_outside_angle_filter'] = \
            #     ids[i].deviation_outside_angle_filter
            # peak_obs.at[i, 'max_m0_abundance'] = \
            #     ids[i].max_m0_abundance
            peak_obs.at[i, 'mzml_path'] = mzml_path
        # Vectorized subtraction
        # peak_obs['envelopes_outside_angle_filter'] = \
        #     peak_obs['envelopes_before_angle_filter'] - \
        #     peak_obs['envelopes_after_angle_filter']

    results = chunk.join(peak_obs)

    return results