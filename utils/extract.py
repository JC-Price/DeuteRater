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

"""
used in the extractor
"""


# TODO: Any temporary values need to be in settings
# TODO: abstract some steps into functions
# TODO: Research how to properly handle exceptions when multiprocessing
# TODO: stress what 'spectrum_index' and 'native_id' mean
# TODO: should we delete invalid envelopes?

# NOTE: this code assumes that an mzml is ordered by retention time
# TODO: should we add logic to check if the mzml is time ordered?


def extract(settings_path, mzml_path, index_to_ID, chunk):
    """Extract data from the mzml according to the identification information
    """
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
    mzml_bounds = dml.get_bounds(mzml_fp, index_to_ID)  # Find what the RT range of the mzML is

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
    hi_rt_bound = chunk.at[len(chunk) - 1, 'rt'] + settings.time_window

    # Search for the scans at the high and low retention time bounds
    lo_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, lo_rt_bound)
    hi_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, hi_rt_bound)
    # if mzml_fp[index_to_ID[hi_spec_idx]] > hi_rt_bound:
    #     hi_spec_idx = hi_spec_idx - 1

    # Logic block for handling out-of-bounds indices
    if lo_spec_idx != -1 and hi_spec_idx != -1:
        # Do nothing if both indices are in bounds
        pass
    elif lo_spec_idx == -1 and hi_spec_idx != -1:
        # If just the higher index is found, assign the lowest index in
        #   the mzml to 'lo_spec_idx'
        lo_spec_idx = mzml_bounds['idx_min']
    elif lo_spec_idx != -1 and hi_spec_idx == -1:
        # If just the lower index is found, assign the highest index in
        #   the mzml to 'hi_spec_idx'
        hi_spec_idx = mzml_bounds['idx_max']
    elif lo_rt_bound < mzml_bounds['rt_min'] < mzml_bounds['rt_max'] < hi_rt_bound:
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

    # TODO: redefine this column as ionMass?
    chunk['mass'] = chunk['mz'] * chunk['z']

    # Instantiate all the identifications in the chunk
    for row in chunk.itertuples(index=True):
        ids.append(
            ID(
                rt=row.rt,
                mz=row.mz,
                mass=row.mass,
                z=row.z,
                n_isos=row.n_isos,
                # cf=row.cf
            )
        )

    # Iterate through all the relevant spectra in the mzml
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
            # catch the exception and move on if the spectrum is not found
            continue

        # only deal with desired ms_level
        if spectrum.ms_level != settings.ms_level:
            continue

        # determine id indices of peaks searches
        # adding and subtracting the floating point error tolerance allows us
        # to include the extremes of the range
        local_window_min = spec_rt - settings.time_window  # + settings.fpe_tolerance)
        local_window_max = spec_rt + settings.time_window  # + settings.fpe_tolerance)
        try:
            lo_slice_index = min(chunk[chunk['rt'] > local_window_min].axes[0].tolist())
            hi_slice_index = max(chunk[chunk['rt'] < local_window_max].axes[0].tolist())
        except:
            continue

        # iterate through relevant ids
        for id in ids[dmt.inclusive_slice(lo_slice_index, hi_slice_index)]:
            charge = id.z
            # instantiate an envelope
            envelope = Envelope(
                peaks=[],
                rt=spec_rt,
                n_lookback=settings.peak_lookback,
                n_lookahead=settings.peak_lookahead
            )

            lo_baseline_lookback = None
            hi_baseline_lookback = None
            lo_baseline_lookahead = None
            hi_baseline_lookahead = None

            peak_range_start = 0 - settings.peak_lookback
            peak_range_end = id.n_isos + settings.peak_lookahead

            # Iterate through all the peaks we want to look for
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
                    lo_baseline_lookback = dmt.find_nearest_index(
                        spec_mzs,
                        id.mz - settings.baseline_lookback
                    )
                    hi_baseline_lookback = index
                    lo_baseline_lookahead = index
                    hi_baseline_lookahead = dmt.find_nearest_index(
                        spec_mzs,
                        id.mz + settings.baseline_lookback
                    )

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
                    envelope.append_peak(Peak(
                        mz=search_mz,
                        # TODO: it might be better to set this to NA
                        abundance=0,
                        i=peak_num
                    ))

            # TODO Do i need to speed this up by removing typecheck?
            # If all of the peaks have been found, add it to the
            #   identification (after determining the baseline)
            # NOTE: baseline is defined as the median abundance of the 100
            #   mz units preceding the m0 peak

            # CQ: Changing baseline to be the MAD of 100 m/z datapoints ahead and behind m0 peak.
            # Adapted from Marginean, I; Tang, K; Smith, RD.; Kelly, R; Picoelectrospray Ionization Mass Spectrometry
            #   Using Narrow-Bore Chemically Etched Emitters, ASMS, 2013

            def mad(values):
                m = median(values)
                return median([abs(a - m) for a in values])

            lookback_baseline = [l for l in spec_abs[dmt.inclusive_slice(lo_baseline_lookback, hi_baseline_lookback)] if
                                 l != 0][-100:]
            lookahead_baseline = [l for l in spec_abs[dmt.inclusive_slice(lo_baseline_lookahead, hi_baseline_lookahead)]
                                  if l != 0][1:101]

            normal_distribution_scale_factor = 1.4826

            # lookback_baseline + lookahead_baseline
            envelope.baseline = normal_distribution_scale_factor * mad(lookback_baseline + lookahead_baseline)

            try:
                id.append_envelope(envelope)
            except Exception as e:
                print("why")
    mzml_fp.close()

    for id in ids:
        id.aggregate_envelopes()

    # TODO: Better variable naming here. obs? I can do better
    # TODO: is there better way to initialize this?
    # TODO: add lookback columns?

    # Initialize the dataframe to send back to the main process
    peak_out = pd.DataFrame(
        index=chunk.index.values,
        columns=['mzs', 'abundances', 'lookback_mzs', 'lookback_abundances',
                 'lookahead_mzs', 'lookahead_abundances', 'rt_min', 'rt_max',
                 'baseline_signal', 'signal_noise', "mads",
                 'mzs_list', 'intensities_list', "rt_list", "baseline_list",
                 'num_scans_combined',
                 'mzml_path']
    )

    # Populate valid rows.
    for row in peak_out.itertuples():
        i = row.Index
        id = ids[i]

        if id.condensed_envelope:
            mzs, abundances = id.condensed_envelope.to_obs()
            lb_mzs, lb_abundances = id.condensed_envelope.lb_obs()
            la_mzs, la_abundances = id.condensed_envelope.la_obs()
            peak_out.at[i, 'mzs'] = mzs
            peak_out.at[i, 'abundances'] = abundances
            peak_out.at[i, 'rt_min'] = id.rt_min
            peak_out.at[i, 'rt_max'] = id.rt_max
            peak_out.at[i, 'baseline_signal'] = id.condensed_envelope.baseline
            peak_out.at[i, 'signal_noise'] = id.signal_noise
            peak_out.at[i, 'lookback_mzs'] = lb_mzs
            peak_out.at[i, 'lookback_abundances'] = lb_abundances
            peak_out.at[i, 'lookahead_mzs'] = la_mzs
            peak_out.at[i, 'lookahead_abundances'] = la_abundances
            peak_out.at[i, 'mads'] = str(id.mads)
            peak_out.at[i, 'num_scans_combined'] = len(id._envelopes)
        if id._unfiltered_envelopes and len([id._unfiltered_envelopes[a] for a in
                                             range(len(id._unfiltered_envelopes))
                                             if id._unfiltered_envelopes[
                                                 a].is_valid]) >= settings.min_envelopes_to_combine:
            mz = [[id._unfiltered_envelopes[k]._peaks[j].mz
                   for j in
                   range(0, len(id._unfiltered_envelopes[k]._peaks))]
                  for k in range(len(id._unfiltered_envelopes))]
            ab = [[id._unfiltered_envelopes[k]._peaks[j].ab
                   for j in
                   range(0, len(id._unfiltered_envelopes[k]._peaks))]
                  for k in range(len(id._unfiltered_envelopes))]
            rt = [id._unfiltered_envelopes[k].rt for k in range(len(id._unfiltered_envelopes))]
            baseline_list = [id._unfiltered_envelopes[k].baseline for k in range(len(id._unfiltered_envelopes))]

            peak_out.at[i, 'mzs_list'] = mz
            peak_out.at[i, 'intensities_list'] = ab
            peak_out.at[i, 'rt_list'] = rt
            peak_out.at[i, 'baseline_list'] = baseline_list

            # Clear the envelopes to save some space. :)
            id._unfiltered_envelopes = None
        peak_out.at[i, 'mzml_path'] = mzml_path

    results = chunk.join(peak_out)

    return results
