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
class and functions for supporting the extractor
"""
import numpy as np
import warnings

try:
    import deuterater.settings as settings
except:
    import DeuteRater.deuterater.settings as settings
from .id import ID


class Molecule(object):
    __slots__ = (
        "name",
        "rt_range",
        "ids",
        "rt_peak_variance",
        "unique_chrom_peaks",
        "reps_with_peak",
        "has_valid",
        "chosen_peak_index",
        "center_scan",
        "error",
        "corresponding_chrom_peaks",
        "corresponding_chrom_peak_max_intensities",
        "why_chosen"
    )

    def __init__(self, name=""):
        self.name = name
        self.ids = dict()
        self.center_scan = None
        self.has_valid = False
        self.rt_range = None
        self.unique_chrom_peaks = list()
        self.chosen_peak_index = None
        self.reps_with_peak = list()
        self.rt_peak_variance = dict()
        self.error = dict()
        self.why_chosen = ""

    def add_id(self, id, rep):
        """
        Adds an id for a specific rep of the given molecule to
        the list of ids available. Determines if there is enough data
        for further peak analysis to be performed.
        """
        self.ids[rep] = id
        if self.center_scan is None:
            self.center_scan = len(id._envelopes) / 2
        if id.is_valid:
            self.has_valid = True

    def analyze_peaks(self, min_window_width=10):
        """
        Looks at chromatographic peaks that appear across all stored ids.
            Finds the peaks that have a unique retention time window and
            stores which reps each peak appears in. The results are stored in the
            unique_chrom_peaks and reps_with_peak attributes.
        
        min_window_width: int
            The minimum number of scans a rt window can be composed of.
        """
        allowed_neutromer_peak_variance = settings.allowed_neutromer_peak_variance

        # Looking at MassHunter across files, it appears a chromatographic peak
        #   has at most 2 minutes of variance in the RT between the most intense
        #   scans. We will use this and use it as a decider on if peaks should
        #   be treated as the same or different peaks across EIC's
        allowed_peak_variance_min = settings.allowed_peak_variance_min
        rt_list = next(iter(self.ids.values()))._get_rt_list()
        whole_rt_range = max(rt_list) - min(rt_list)
        distance_between_scans = whole_rt_range / len(rt_list)
        allowed_scan_variance = int(allowed_peak_variance_min / distance_between_scans)

        # Grab peaks from all the ids for comparison. Also, remove any peaks
        #   that are at the beginning or ending scan of the 3-minute window.
        from copy import deepcopy, copy
        available_chrom_peaks = list()
        id_neutromer_peak_indexes = list()
        chrom_peak_indexes = list()
        rep_list = list()
        for rep, id in self.ids.items():
            # id.divide_chromatography(True)
            rep_list.append(rep)
            id_neutromer_peak_indexes.append(id.neutromer_peak_maximums)
            available_chrom_peaks.append(id.rt_windows)
            chrom_peak_indexes.append(id.rt_peak_index)
            available_chrom_peaks[-1].insert(0, rep)
            chrom_peak_indexes[-1].insert(0, rep)

            # Remove any peaks that has its center at the end or beginning of 3-minute window.
            if len(id.rt_peak_index) > 1:
                end_range = None
                start_range = None
                if id.rt_peak_index[1] == 0:
                    start_range = available_chrom_peaks[-1].pop(1)
                    chrom_peak_indexes[-1].pop(1)
                if id.rt_peak_index[-1] == (len(id._envelopes) - 1):
                    end_range = available_chrom_peaks[-1].pop(-1)
                    chrom_peak_indexes[-1].pop(-1)
                for m in id_neutromer_peak_indexes[-1]:
                    if start_range and m and m[0] in range(start_range[0], start_range[1] + 1):
                        m.pop(0)
                    if end_range and m and m[-1] in range(end_range[0], end_range[1] + 1):
                        m.pop(-1)

        # Remove peaks that are not wide enough (#Scans < 10)
        #   Also remove any windows that are composed by only
        #   1 or 2 neutromers.
        for id_num in range(len(available_chrom_peaks)):
            peak_windows = deepcopy(available_chrom_peaks[id_num])
            peak_maximums = deepcopy(chrom_peak_indexes[id_num])
            neutromer_peaks = id_neutromer_peak_indexes[id_num]
            for peak_num in range(1, len(peak_windows)):
                try:
                    peak_window = peak_windows[peak_num]
                    peak_max = peak_maximums[peak_num]
                    id = self.ids[rep_list[id_num]]
                except:
                    continue

                # For a peak to be considered, it must be comprised of all 3 neutromer peaks, and the
                #   peak of each neutromer's elution curve must be within +-5 scans of the combined peak.
                all_neutromers_exist = all(
                    [any(
                        [peak_window[0] <= neutromer_peaks[neutromer_num][p] <= peak_window[1] and
                         peak_max - allowed_neutromer_peak_variance <= neutromer_peaks[neutromer_num][p]
                         <= peak_max + allowed_neutromer_peak_variance
                         for p in range(len(neutromer_peaks[neutromer_num]))]
                    ) for neutromer_num in range(len(neutromer_peaks))])
                num_valid_in_window = len(
                    [envelope.is_valid for envelope in id._envelopes[peak_window[0]:peak_window[1]] if
                     envelope.is_valid])
                if num_valid_in_window < min_window_width or not all_neutromers_exist:
                    available_chrom_peaks[id_num].remove(peak_windows[peak_num])
                    chrom_peak_indexes[id_num].remove(peak_maximums[peak_num])

        # Compare peaks across reps to find the number of unique peaks.
        index_checking_list = copy(chrom_peak_indexes)
        index_checking_list = [peak for peak in index_checking_list if len(peak) > 1]
        window_checking_list = copy(available_chrom_peaks)
        window_checking_list = [peak_list for peak_list in window_checking_list if len(peak_list) > 1]

        left_positions = list()
        right_positions = list()
        where_peak_appears = list()
        peaks_that_appear = list()
        while len(window_checking_list) != 0:
            # Add the peak with the smallest starting RT:
            index_positions = list()
            current_layer_windows = [peak_list[1] for peak_list in window_checking_list]
            current_layer_peak_indexes = [peak_list[1] for peak_list in index_checking_list]
            index_to_add = np.argmin([layer_peak[0] for layer_peak in current_layer_windows])
            left_positions.append(current_layer_windows[index_to_add][0])
            right_positions.append(current_layer_windows[index_to_add][1])
            index_positions.append(current_layer_peak_indexes[index_to_add])
            where_peak_appears.append(list())
            peaks_that_appear.append(dict())
            # Expand the window and count how many reps have this peak:
            for list_index in range(len(window_checking_list)):
                current_window = window_checking_list[list_index][1]
                current_index = index_checking_list[list_index][1]
                # If the peak has a starting position that is larger than current end of the peak,
                #   it is a new one and should be skipped.
                if right_positions[-1] < current_window[0] or (
                        min(index_positions) + allowed_scan_variance) < current_index:
                    continue

                # If the end of the peak is at a later RT but has a starting RT within the current window,
                #   just expand the RT end.
                elif left_positions[-1] < current_window[0] and right_positions[-1] < current_window[1]:
                    right_positions[-1] = current_window[1]

                # Add the current rep to the list of reps the peak appears in
                where_peak_appears[-1].append(window_checking_list[list_index][0])
                peaks_that_appear[-1][window_checking_list[list_index][0]] = current_window
                index_positions.append(current_index)
                window_checking_list[list_index].remove(current_window)
                index_checking_list[list_index].remove(current_index)

            # Remove all empty lists:
            window_checking_list = [peak_list for peak_list in window_checking_list if len(peak_list) > 1]
            index_checking_list = [peak_list for peak_list in index_checking_list if len(peak_list) > 1]

        unique_peaks = [(left_positions[i], right_positions[i]) for i in range(len(left_positions))]
        self.reps_with_peak = where_peak_appears
        self.unique_chrom_peaks = unique_peaks
        self.corresponding_chrom_peaks = peaks_that_appear

        # TODO: figure out if this would be useful to include - Ben D
        # Remove peaks that switch neutromer peak order when
        #   sorted by intensity.
        index_to_remove = []
        for num in range(len(self.corresponding_chrom_peaks)):
            peak = self.corresponding_chrom_peaks[num]
            if len(peak) < 2:
                continue
            # neutromer_intensity_lists = [[[str(ab.ab) for ab in envelope._peaks[1:4]] for envelope in self.ids[addu]._envelopes[rt_range[0]: rt_range[1] + 1]] for addu, rt_range in peak.items()]
            # unique_orders = np.unique([intensities for intensities in neutromer_intensity_lists])
            unique_order = [np.unique([",".join(np.argsort([ab.ab for ab in envelope._peaks[1:4]]).astype(str))
                                       for envelope in self.ids[addu]._envelopes[int((rt_range[1]+rt_range[0])/2 - 3): int((rt_range[1]+rt_range[0])/2 + 4)]]
                                      , return_counts=True) for addu, rt_range in peak.items()]
            unique_orders = np.unique([i[0][np.argmax(i[1])] for i in unique_order])
            if len(unique_orders) > 1:
                index_to_remove.insert(0, num)

        for index in index_to_remove:
            self.reps_with_peak.pop(index)
            self.unique_chrom_peaks.pop(index)
            self.corresponding_chrom_peaks.pop(index)

        self._find_peak_variance()

        # Store max height of each chromatographic peak.
        stored_intensities = list()
        for peak in self.corresponding_chrom_peaks:
            curr_id_list = dict()
            for rep, rt_range in peak.items():
                id = self.ids[rep]
                intensities = np.sum(np.array(id._get_peak_list()), axis=0)
                curr_id_list[rep] = max(intensities[rt_range[0]:rt_range[1]])
            stored_intensities.append(curr_id_list)
        self.corresponding_chrom_peak_max_intensities = stored_intensities

    def _find_peak_variance(self):
        rep_variance = dict()
        for rep, id in self.ids.items():
            # We can only look at neutromer peak variance if the id is valid.
            if not id.is_valid or len(self.reps_with_peak) == 0:
                continue
            # Used to plot the current ID's chromatography. Remove after testing

            unique_peaks = self.unique_chrom_peaks
            neutromer_peaks = id.neutromer_peak_maximums

            peak_variance = list()

            total_peaks_per_neutromer = [len(neutromer_peaks[i]) - 1 for i in range(id.n_isos)]
            # If any neutromers have no peaks, then no peaks will
            # work for this given rep
            if -1 in total_peaks_per_neutromer:
                for peak in unique_peaks:
                    peak_variance.append([np.nan])
                rep_variance[rep] = peak_variance
                continue

            current_neutromer_peaks_index = [0 for i in range(id.n_isos)]
            current_unique_peak_index = 0

            more_windows_to_compare = True
            while more_windows_to_compare:
                # Continue cycling though peaks available for a given rep until one appears
                if rep not in self.reps_with_peak[current_unique_peak_index]:
                    peak_variance.append([np.nan])
                    current_unique_peak_index += 1
                    if current_unique_peak_index == len(unique_peaks):
                        break
                    continue
                valid_window = True

                # Get information to analyze. This will change every iteration of the while loop
                current_neutromer_peaks = [neutromer_peaks[i][current_neutromer_peaks_index[i]] for i in
                                           range(id.n_isos)]
                current_unique_peak = unique_peaks[current_unique_peak_index]
                # Check if all neutromer peak maximums are in
                #   the current peak's window.
                for i in range(id.n_isos):
                    current_neutromer_peak = current_neutromer_peaks[i]
                    if current_unique_peak[0] < current_neutromer_peak < current_unique_peak[1]:
                        # The neutromer peak is inside the window.
                        continue
                    elif current_unique_peak[0] > current_neutromer_peak:
                        # The current window is later than the current neutromer peak,
                        #   so I need to move it forward one index, if possible. If it
                        #   is the last one, then stop the variance analysis, because
                        #   nothing else can be done.
                        valid_window = False
                        if current_neutromer_peaks_index[i] == total_peaks_per_neutromer[i]:
                            more_windows_to_compare = False
                            break
                        else:
                            current_neutromer_peaks_index[i] += 1
                    elif current_unique_peak[1] < current_neutromer_peak:
                        # The neutromer_peak is later than the current peak, so this
                        #   unique peak should be skipped. If it is the last unique peak,
                        #   then neutromer peak variance calculation should stop.
                        valid_window = False

                        # Add None to show that that peak does not appear for all
                        #   neutromers and thus should not be treated as valid.
                        # TODO: Consider removing the rep being analyzed from the file the given peak appears in.
                        # self.reps_with_peak[current_unique_peak_index].remove(rep)
                        peak_variance.append([np.nan])

                        if current_unique_peak_index == len(self.unique_chrom_peaks) - 1:
                            more_windows_to_compare = False
                            break
                        else:
                            current_unique_peak_index += 1
                            break
                if valid_window:
                    # neutromer_peak_variance = max(current_neutromer_peaks) - min(current_neutromer_peaks)
                    neutromer_peak_variance = current_neutromer_peaks
                    peak_variance.append(neutromer_peak_variance)
                    current_unique_peak_index = current_unique_peak_index + 1
                    current_neutromer_peaks_index = [current_neutromer_peaks_index[i] + 1 for i in
                                                     range(len(current_neutromer_peaks_index))]
                    if current_unique_peak_index == len(unique_peaks):
                        break
                    if any([current_neutromer_peaks_index[i] == total_peaks_per_neutromer[i] + 1 for i in
                            range(id.n_isos)]):
                        break

            # Make note that future peaks don't occur in the current rep.
            while len(peak_variance) != len(self.unique_chrom_peaks):
                peak_variance.append([np.nan])
            rep_variance[rep] = peak_variance
        self.rt_peak_variance = rep_variance

    def choose_peak(self):
        if len(self.unique_chrom_peaks) == 0:
            return
        # Grab the earliest peak (default)
        chosen_peak = 0

        num_reps_with_peak = np.array([-len(peak) for peak in self.reps_with_peak])
        most_peaks = max([len(peak) for peak in self.reps_with_peak])
        num_reps_percentage = np.array([len(peak) / most_peaks for peak in self.reps_with_peak])
        with warnings.catch_warnings():
            # nanmean will return a RuntimeWarning if I am calculating a mean on all np.nan.
            #   There is no other reason for a RuntimeWarning here, so I will suppress it
            #   because the warning is telling me something I already know and am handling correctly.
            warnings.simplefilter("ignore", category=RuntimeWarning)

            def stat_range(a):
                return max(a) - min(a)

            def stat_range_percent(a):
                return 1 - (max(a) - min(a)) / settings.allowed_neutromer_peak_variance

            variance_totals = np.array([np.nanmean(
                [stat_range(rep[peak]) for rep in self.rt_peak_variance.values() if not all(np.isnan(rep[peak]))]) for
                                        peak in range(len(self.unique_chrom_peaks))])
            variance_percentages = np.array([np.nanmean(
                [stat_range_percent(rep[peak]) for rep in self.rt_peak_variance.values() if
                 not all(np.isnan(rep[peak]))]) for peak in range(len(self.unique_chrom_peaks))])
        dist_from_center = np.array(
            [abs(((peak[0] + peak[1]) / 2) - self.center_scan) for peak in self.unique_chrom_peaks])
        dist_from_center_percentage = np.array(
            [1 - abs(((peak[0] + peak[1]) / 2) - self.center_scan) / self.center_scan for peak in
             self.unique_chrom_peaks])

        stored_intensities = self.corresponding_chrom_peak_max_intensities
        combined_intensities = [sum(peak.values()) for peak in stored_intensities]
        intensity_percentage = np.array([number / max(combined_intensities) for number in combined_intensities])

        try:
            combined_percentages = settings.adduct_weight * num_reps_percentage + \
                                   settings.variance_weight * variance_percentages + \
                                   settings.ID_weight * dist_from_center_percentage + \
                                   settings.intensity_weight * intensity_percentage
        except:
            print("percentage error")
        import scipy.stats as ss
        # rep_rank = ss.rankdata(num_reps_with_peak, method="dense")
        # variance_rank = ss.rankdata(variance_totals, method="dense")
        # rt_rank = ss.rankdata(dist_from_center, method="dense")
        # combined_rank = rep_rank + variance_rank + rt_rank
        # best_peak = np.argwhere(combined_rank == min(combined_rank))[0]
        best_percentage = np.argwhere(combined_percentages == max(combined_percentages))[0]
        if len(best_percentage) != 1:
            chosen_peak = np.argwhere(dist_from_center == min(dist_from_center))[0][0]
            self.why_chosen += "TIE: RT Used\n"
        else:
            chosen_peak = best_percentage[0]
        self.why_chosen += f"{round(num_reps_percentage[chosen_peak], 4)}, " \
                           f"{round(variance_percentages[chosen_peak], 4)}, " \
                           f"{round(dist_from_center_percentage[chosen_peak], 4)}, " \
                           f"{round(intensity_percentage[chosen_peak], 4)}"

        # if len(np.argwhere(num_reps_with_peak == max(num_reps_with_peak))) == 1:
        #     chosen_peak = np.argwhere(num_reps_with_peak == max(num_reps_with_peak))[0][0]
        #     self.why_chosen = "number of reps"
        # elif len(np.argwhere(variance_totals == min(variance_totals))) == 1:
        #     chosen_peak = np.argwhere(variance_totals == min(variance_totals))[0][0]
        #     self.why_chosen = "np variance"
        # else:
        #     chosen_peak = np.argwhere(dist_from_center == min(dist_from_center))[0][0]
        #     self.why_chosen = "rt position"
        try:
            self.rt_range = self.unique_chrom_peaks[chosen_peak]
        except:
            print("range error")
        self.chosen_peak_index = chosen_peak

    def extract_with_chosen_peak(self):
        if self.chosen_peak_index == None:
            for rep, id in self.ids.items():
                self.error[rep] = "No peak was able to be chosen"
            return
        for rep, id in self.ids.items():
            # id.divide_chromatography(True)
            if self.chosen_peak_index is not None and rep in self.reps_with_peak[self.chosen_peak_index]:
                rt_range = self.corresponding_chrom_peaks[self.chosen_peak_index][rep]
                id._envelopes = id._envelopes[rt_range[0]: rt_range[1] + 1]
                id.aggregate_envelopes()
                if not id.condensed_envelope:
                    self.error[rep] = "CFS, Angle, Threshold Filtered resulted in less than 10 scans."
            else:
                self.error[rep] = "Chromatographic Peak does not exist for this rep."

    def update_output_file(self, df):
        # df_w_index = df.reset_index(drop=True)
        df["row_num"] = np.arange(0, df.shape[0])
        df.set_index("row_num", drop=False, inplace=True)
        for rep, id in self.ids.items():
            name_column = list(set(df.columns).intersection(set(["Lipid Name", "Sequence"])))
            row = df[df[name_column[0]] == self.name]
            row = row.loc[row["name_check"] == rep]
            try:
                row_index = row.iloc[0]['row_num']
            except:
                print("There is an errror")
                continue
            if id.condensed_envelope:
                mzs, abundances = id.condensed_envelope.to_obs()
                lb_mzs, lb_abundances = id.condensed_envelope.lb_obs()
                la_mzs, la_abundances = id.condensed_envelope.la_obs()
                df.at[row_index, 'mzs'] = mzs
                df.at[row_index, 'abundances'] = abundances
                df.at[row_index, 'rt_min'] = id.rt_min
                df.at[row_index, 'rt_max'] = id.rt_max
                df.at[row_index, "Extraction_Updated"] = self.why_chosen
                df.at[row_index, 'baseline_signal'] = id.condensed_envelope.baseline
                # df.at[row_index, 'signal_noise'] = ""
                df.at[row_index, 'signal_noise'] = str(id.signal_noise)
                df.at[row_index, 'lookback_mzs'] = lb_mzs
                df.at[row_index, 'lookback_abundances'] = lb_abundances
                df.at[row_index, 'lookahead_mzs'] = la_mzs
                df.at[row_index, 'lookahead_abundances'] = la_abundances
                df.at[row_index, 'mads'] = str(id.mads)
                df.at[row_index, 'num_scans_combined'] = len(id._envelopes)
            else:
                try:
                    df.at[row_index, 'mzs'] = np.nan
                except:
                    print("There is an errror in the else")
                    continue
                df.at[row_index, 'abundances'] = np.nan
                df.at[row_index, 'rt_min'] = np.nan
                df.at[row_index, 'rt_max'] = np.nan
                df.at[row_index, 'baseline_signal'] = np.nan
                df.at[row_index, 'signal_noise'] = np.nan
                df.at[row_index, 'lookback_mzs'] = np.nan
                df.at[row_index, 'lookback_abundances'] = np.nan
                df.at[row_index, 'lookahead_mzs'] = np.nan
                df.at[row_index, 'lookahead_abundances'] = np.nan
                df.at[row_index, 'mads'] = np.nan
                df.at[row_index, 'num_scans_combined'] = 0
                df.at[row_index, "Extraction_Updated"] = "no peak used"
                df.at[row_index, "Extraction_Error"] = self.error[rep]
        return df
