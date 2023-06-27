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
import sys
from pymzml import pymzml
import warnings
import pandas as pd
from numpy import median

import time

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings
import utils.mzml as dml
import utils.mathDR as dmt
from utils.exc import EmptyIdChunkWarning
from constants import NEUTRON

from obs.peak import Peak
from obs.envelope import Envelope
from obs.id import ID


# TODO: Any temporary values need to be in settings
# TODO: abstract some steps into functions
# TODO: Research how to properly handle exceptions when multiprocessing
# TODO: stress what 'spectrum_index' and 'native_id' mean
# TODO: should we delete invalid envelopes?

# NOTE: this code assumes that an mzml is ordered by retention time
# TODO: should we add logic to check if the mzml is time ordered?


def extract(settings_path, mzml_path, id_path, is_ims=True, index_to_ID=None, chunk=None):
    # Extract data from the mzml according to the identification information

    # A rough outline of how the logic flows.
    # EXTRACT:
    # for each scan in the reader
    #   in ms level 1
    #   for each id in the window
    #       isotope extraction specific logic

    # first=True

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

    if chunk is None:
        chunk = pd.read_csv(id_path)  # this comes from the csv
        # print('chunk ==', chunk) #remove this after debugging
        chunk = chunk.sort_values(by=['rt_minutes', 'dt_msec'], ignore_index=True)

    # BD: Don't quite understand what the purpose of index_to_ID is?
    if index_to_ID is None:
        print("starting to index")
        start = time.time()
        index_to_ID = {i: i for i in range(2918159)}  # D0
        index_to_ID = {i: i for i in range(2917079)}  # D32
        # index_to_ID = {spec.index: spec.ID for spec in mzml_fp}
        print("Time to index: {}".format(time.time() - start))
    mzml_bounds = dml.get_bounds(mzml_fp, index_to_ID)

    # mzml_bounds = {
    # 	'idx_min': 0,
    # 	'idx_max': mzml_fp.get_spectrum_count() - 1,
    # 	'rt_min': mzml_fp[0].scan_time_in_minutes(),
    # 	'rt_max': mzml_fp[mzml_fp.get_spectrum_count() - 1].scan_time_in_minutes()
    # }

    # index_to_ID = {x: x for x in range(mzml_bounds["idx_max"])}

    # mzml_bounds = dml.get_bounds(mzml_fp, index_to_ID)  # Find what the RT range of the mzML is

    # Check for an empty id file chunk
    if len(chunk) <= 0:
        warnings.warn(EmptyIdChunkWarning(
            'There are no identifications in this chunk'
        ))

    # Set the high and low retention time bounds, based on the chunk of the
    #   identification file
    lo_rt_bound = chunk.loc[0, 'rt_minutes'] - settings.time_window
    if lo_rt_bound < 0:
        lo_rt_bound = 0
    hi_rt_bound = chunk.loc[len(chunk) - 1, 'rt_minutes'] + settings.time_window
    # print('lower bound==', lo_rt_bound)  # remove after debugging
    # print('upper bound==', hi_rt_bound)  # remove after debugging
    # print('mzml_lower==', mzml_fp[min(index_to_ID)].scan_time_in_minutes())  # remove after debugging
    # print('mzml_higher==', mzml_fp[max(index_to_ID)].scan_time_in_minutes())  # remove after debugging
    # Search for the scans at the high and low retention time bounds
    lo_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, lo_rt_bound)  # remove hags1 and hags2 after debugging
    hi_spec_idx = dml.retention_time_search(mzml_fp, index_to_ID, hi_rt_bound)

    #    if mzml_fp[index_to_ID[hi_spec_idx]] > hi_rt_bound:
    #        hi_spec_idx = hi_spec_idx - 1
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

    # Instantiate all the identifications in the chunk
    for row in chunk.itertuples(index=True, name="ID"):
        new_id = ID(
            rt=row.rt_minutes,
            mz=row.mz,
            mass=row.mass,
            z=row.z,
            n_isos=row.n_isos,
        )
        if is_ims:
            new_id.dt = row.dt_msec
        ids.append(new_id)

    # Iterate through all the relevent spectrums in the mzml
    for spectrum_index in dmt.inclusive_range(lo_spec_idx, hi_spec_idx):
        # apply the index_to_ID map in order to access the correct spectrum
        native_id = index_to_ID[spectrum_index]
        try:
            # try to access this spectrum
            spectrum = mzml_fp[native_id]
            spec_rt = spectrum.scan_time_in_minutes()
            if is_ims:
                spec_dt = dml.get_dt(spectrum)
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
        local_window_min = \
            spec_rt - settings.time_window  # + settings.fpe_tolerance)
        local_window_max = \
            spec_rt + settings.time_window  # + settings.fpe_tolerance)

        if is_ims:
            drift_window = 0.5
            local_drift_min = spec_dt - drift_window
            local_drift_max = spec_dt + drift_window
            if 33.5 > spec_dt > 33.4 and 36.1 < spec_rt < 36.3:
                print('hi')
            try:
                lo_slice_index = \
                    min(chunk[(chunk['rt_minutes'] > local_window_min) & (chunk['dt_msec'] > local_drift_min)].axes[
                            0].tolist())
                hi_slice_index = \
                    max(chunk[(chunk['rt_minutes'] < local_window_max) & (chunk['dt_msec'] < local_drift_max)].axes[
                            0].tolist())
            except Exception as e:
                continue
        else:
            try:
                lo_slice_index = \
                    min(chunk[(chunk['rt_minutes'] > local_window_min)].axes[0].tolist())
                hi_slice_index = \
                    max(chunk[(chunk['rt_minutes'] < local_window_max)].axes[0].tolist())
            except Exception:
                continue

        # iterate through relevant ids
        for id in ids[dmt.inclusive_slice(lo_slice_index, hi_slice_index)]:
            charge = id.z
            # if is_ims and id.dt < spec_dt - drift_window or id.dt > spec_dt + drift_window:
            # if first:
            # 	print("FIX DT ID RANGE CHECKING!!!")
            # 	first = False
            # continue
            # instatntiate an envelope
            envelope = Envelope(
                peaks=[],
                rt=spec_rt,
                n_lookback=settings.peak_lookback,
                n_lookahead=settings.peak_lookahead
            )
            if is_ims:
                envelope.dt = spec_dt

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
                if len(spec_mzs) > 1 and abs(spec_mzs[index] - search_mz) < reach:
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
            envelope.baseline = normal_distribution_scale_factor * mad(
                lookback_baseline + lookahead_baseline) * 3  # lookback_baseline + lookahead_baseline  #

            try:
                id.append_envelope(envelope)
            except Exception as e:
                print("why")
    mzml_fp.close()

    for id in ids:
        id.aggregate_envelopes()

    # TODO: Better variable naming here. obs? I can do better
    # TODO: is there better way to initialize this?
    # cols for peak_out(1) ('mzs', 'abundances', 'm-1_mz', 'm-1_abundance','m_end+1_mz', 'm_end+1_abundance', 'rt_min', 'rt_max', 'baseline_signal', 'signal_2_noise', "mads", "num_scans_combined")
    cols = ['mzs', 'abundances', 'm-1_mz', 'm-1_abundance',
            'm_end+1_mz', 'm_end+1_abundance', 'rt_min', 'rt_max',
            'baseline_signal', 'signal_2_noise', "mads",
            'num_scans_combined', 'summed_signal']
    # cols for peak_out (2)('mzs_list', 'intensities_list', "rt_list", "baseline_list", 'id_path', 'mzml_path', 'summed_signal')
    cols_2 = ['mzs_list', 'intensities_list', "rt_list", "baseline_list", 'id_path', 'mzml_path']

    if is_ims:
        cols_2.append("dt_list")
    # Initialize the dataframe to send back to the main process
    peak_out = pd.DataFrame(
        index=chunk.index.values,
        columns=cols
    )
    peak_out_2 = pd.DataFrame(
        index=chunk.index.values,
        columns=cols_2)
    # Populate valid rows.
    for row in peak_out.itertuples():
        i = row.Index
        id = ids[i]
        # print(ids)
        if id.condensed_envelope:
            mzs, abundances = id.condensed_envelope.to_obs()
            summed_signal = str(id.total_signal)
            signal_2_noise = [float(a) / float(id.condensed_envelope.baseline) for a in abundances]
            lb_mzs, lb_abundances = id.condensed_envelope.lb_obs()
            la_mzs, la_abundances = id.condensed_envelope.la_obs()
            peak_out.loc[i, 'mzs'] = mzs
            peak_out.loc[i, 'abundances'] = abundances
            peak_out.loc[i, 'summed_signal'] = summed_signal
            peak_out.loc[i, 'rt_min'] = id.rt_min
            peak_out.loc[i, 'rt_max'] = id.rt_max
            peak_out.loc[i, 'baseline_signal'] = id.condensed_envelope.baseline
            peak_out.loc[i, 'signal_2_noise'] = signal_2_noise
            peak_out.loc[i, 'm-1_mz'] = lb_mzs
            peak_out.loc[i, 'm-1_abundance'] = lb_abundances
            peak_out.loc[i, 'm_end+1_mz'] = la_mzs
            peak_out.loc[i, 'm_end+1_abundance'] = la_abundances
            peak_out.loc[i, 'mads'] = str(id.mads)
            peak_out.loc[i, 'num_scans_combined'] = len(id.envelopes)
        # if is_ims:
        #     peak_out.loc[i, 'dt_min'] = id.dt_min
        #     peak_out.loc[i, 'dt_max'] = id.dt_min

        # BD: Had some issues with PyCharm type checking... Added # noqa to ignore warnings for id.unfiltered_envelopes
        # Also suppressed unresolved reference warnings for the entire statement
        # noinspection PyUnresolvedReferences
        if id.unfiltered_envelopes and len([id.unfiltered_envelopes[a] for a in
                                            range(len(id.unfiltered_envelopes)) # noqa
                                            if id.unfiltered_envelopes[
                                                a].is_valid]) >= settings.min_envelopes_to_combine:
            mz = [[id.unfiltered_envelopes[k]._peaks[j].mz  # noqa
                   # where mz is defined, a long list of lists. #This seems to have worked well. Try taking away the outer brackets
                   for j in
                   range(0, len(id.unfiltered_envelopes[k]._peaks))]  # noqa
                  for k in range(len(id.unfiltered_envelopes))] # noqa
            ab = [[id.unfiltered_envelopes[k]._peaks[j].ab  # noqa # this one did not seem to work
                   for j in
                   range(0, len(id.unfiltered_envelopes[k]._peaks))]  # noqa
                  for k in range(len(id.unfiltered_envelopes))] # noqa
            rt = [id.unfiltered_envelopes[k].rt for k in range(len(id.unfiltered_envelopes))] # noqa
            baseline_list = [id.unfiltered_envelopes[k].baseline for k in range(len(id.unfiltered_envelopes))] # noqa

            if is_ims:
                dt = [id.unfiltered_envelopes[k].dt for k in range(len(id.unfiltered_envelopes))] # noqa
                peak_out_2.loc[i, 'dt_list'] = dt  # [1,2,3,4,5] #change back to dt after debugging
            peak_out_2.loc[i, 'mzs_list'] = mz  # change back to mz after debugging
            peak_out_2.loc[i, 'intensities_list'] = ab  # [1,2,3,4,5] #change back to ab after debugging
            peak_out_2.loc[i, 'rt_list'] = rt  # [1,2,3,4,5] #change back to rt after debugging
            peak_out_2.loc[
                i, 'baseline_list'] = baseline_list  # [1,2,3,4,5] #change back to baseline_list after debugging
            print(mz)
            print(type(mz))

            # Clear the envelopes to save some space. :)
            id.unfiltered_envelopes = None
        peak_out_2.loc[i, 'mzml_path'] = mzml_path
        peak_out_2.loc[i, 'id_path'] = id_path
    peak_out_f = peak_out.join(peak_out_2)
    results = chunk.join(peak_out_f)
    results.to_csv('results_test_test.csv', sep="\t", index=False)  # remove after debugging
    results['id_path'] = results["id_path"].apply(lambda x: x.replace('\\t', '/t'))
    results['mzml_path'] = results["mzml_path"].apply(lambda x: x.replace('\\t', '/t'))
    print(results)  # remove after debugging  #this works fine,why is it not able to print after returning?

    return results, chunk, peak_out_f, peak_out_2, peak_out  # do not need to return chunk and peakout(s)


# ToDo: Fix main
# def main():
# 	id_path = 'D:/IMS_ID_File_PE36_2.csv'
# 	mzml_path = "D:/A3D32M1_pos_3.20-end_DI3.d.HRdm.mzML"
# 	df = extract("settings.yaml", mzml_path, id_path, is_ims=True)
# 	# df['id_path'] = df["id_path"].apply(lambda x: x.replace('\\t', '/t'))
# 	# df['mzml_path'] = df["mzml_path"].apply(lambda x: x.replace('\\t', '/t'))
# 	# print(df) # remove this when finished debugging
#
# 	if df != -1:
# 		print(df)
# 		df.to_csv("test1_only_extracted_sorted32_5dtw.csv", sep="\t")
# 	else:
# 		raise 'Failed to generate exportable dataframe because requirements for "if or elif" statements were not met in extract function'


# if __name__ == "__main__":
# main()

def main():
    paths = sys.argv[1:]
    # mzml_path = '/home/JCPriceLab/Desktop/DeuteRater Testing/A4_PQC_3_pos_3.20-end_DI3.d.DeMP.csv'
    # id_path = '/home/JCPriceLab/Desktop/DeuteRater Testing/IMS_ID_File_PE36_2.csv'
    # results, chunk, peak_out_f, peak_out_2, peak_out = extract("settings.yaml", mzml_path, id_path, is_ims=True)

    mzml_path = paths[0]
    id_path = paths[1]

    results, chunk, peak_out_f, peak_out_2, peak_out = extract("settings.yaml", mzml_path, id_path, is_ims=True)

# chunk.to_csv('chunk_test_output.tsv', sep="\t", index=False)
# peak_out_f.to_csv('peak_out_f_output.tsv', sep="\t", index=False)
# peak_out_2.to_csv('peak_out_2_output.tsv', sep='\t', index=False)
# peak_out.to_csv('peak_out_1_output.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()
