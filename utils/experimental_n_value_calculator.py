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

try:
    from utils.n_value_calc_emass import emass
except ModuleNotFoundError as e:
    from n_value_calc_emass import emass

from functools import partial
import pandas as pd
import multiprocessing as mp
import concurrent.futures as cf
from tqdm import tqdm
import deuterater.settings as settings
import numpy as np
import traceback
import re
import matplotlib.pyplot as plt

output_columns = ['n_value', 'n_value_stddev', "num_nv_time_points", "cv", 'Median_FN_stddev', 'n_val_lower_margin',
                  'n_val_upper_margin', 'All_n_values', 'N_value_time_points', 'Filtered_out_N_values',
                  'Filtered_out_N_value_time_points', 'Reason_for_filtering']
error_statement = ["error occurred", "error occurred", "error occurred", "error occurred", "error occurred",
                   "error occurred", "error occurred", "error occurred", "error occurred", "error occurred", "error occurred",
                   "error occurred"]
timepoints_statement = ["no valid time points", 'n_value_stddev', "num_nv_time_points", "cv", 'Median_FN_stddev',
                        'n_val_lower_margin', 'n_val_upper_margin', 'All_n_values', 'N_value_time_points',
                        'Filtered_out_N_values', 'Filtered_out_N_value_time_points', 'Reason_for_filtering']


def bootstrap_median_ci(data, n_iterations=1000, ci=95, seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Storage for bootstrapped medians
    medians = np.empty(n_iterations)

    # Bootstrap resampling
    for i in range(n_iterations):
        resample = np.random.choice(data, size=len(data), replace=True)
        medians[i] = np.median(resample)

    # Compute confidence interval bounds
    lower_percentile = (100 - ci) / 2
    upper_percentile = 100 - lower_percentile

    lower_bound = np.percentile(medians, lower_percentile)
    upper_bound = np.percentile(medians, upper_percentile)

    # Compute the median of the original dataset
    median_value = np.median(data)

    # Compute the asymmetric margins
    lower_margin = median_value - lower_bound
    upper_margin = upper_bound - median_value

    return lower_margin, upper_margin


class NValueCalculator:
    def __init__(self, dataframe, settings_path, biomolecule_type, out_path, graphs_location=None):
        # this DataFrame should only contain the correct columns
        settings.load(settings_path)
        self.settings_path = settings_path
        self.biomolecule_type = biomolecule_type
        self.out_path = out_path
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.75)
        else:
            self._n_processors = settings.n_processors
        # breaks windows/python interactions if too many cores are used.  very niche application but still relevant
        if self._n_processors > 60:
            self._n_processors = 60

        self.full_df = dataframe
        self.output_columns = output_columns
        self.prepared_df = None

        if settings.graph_n_value_calculations:
            self.graph_folder = graphs_location
        else:
            self.graph_folder = None

    def prepare_model(self):
        if "Adduct_cf" in self.full_df.columns:
            df = self.full_df.rename(columns={
                'Adduct_cf': 'adduct_chemical_formula',
                'cf': 'chemical_formula',
                'abundances': 'intensities',
                'labeling': 'enrichment'})
        else:
            # Renames columns to proper names depending on if a lipid or protein is being passed in
            df = self.full_df.rename(columns={
                'cf': 'chemical_formula',
                'abundances': 'intensities',
                'labeling': 'enrichment'})
            df["adduct_chemical_formula"] = df["chemical_formula"]

        df['intensities'] = df['intensities'].astype(str).apply(lambda x: x[1:-1])
        df['chemical_formula'] = df['chemical_formula'].astype(str)

        # Normalize intensities
        for row in df.itertuples():  # Normalise the intensity list by sum
            temp = [float(x) for x in row.intensities.split(', ')]
            s = sum(temp)
            df.at[row.Index, 'intensities'] = [float(i) / s for i in temp]

        df = df.reset_index(drop=False)

        if self.biomolecule_type == "Peptide":
            df = df[['index', 'Protein ID', 'chemical_formula', 'intensities', 'enrichment', 'sample_group',
                     'bio_rep', 'calculate_n_value', 'time']]
            df.loc[:, "divider"] = df['Protein ID'] + df['sample_group']
        else:
            df = df[['index', 'Lipid Unique Identifier', 'chemical_formula', 'adduct_chemical_formula',
                     'intensities', 'enrichment', 'sample_group', 'bio_rep', 'calculate_n_value', 'time']]
            df.loc[:, "divider"] = df['Lipid Unique Identifier'] + df['sample_group']

        df.sort_values(by="divider")

        self.prepared_df = df

    def run(self):
        self.prepare_model()

        groups = self.prepared_df.groupby('divider', sort=False)

        nv_function = partial(NValueCalculator.analyze_group, self, save_data=settings.save_n_value_data,
                              make_graphs=settings.graph_n_value_calculations)

        results = list()
        n_value_data = list()

        # if either of these settings are turned on, we need to run this step through debug mode so we can get all the results back
        if settings.graph_n_value_calculations or settings.save_n_value_data:
            settings.debug_level = 1

        # settings.debug_level = 1
        n_data = []
        if settings.debug_level == 0:
            with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                results = list(
                    tqdm(executor.map(nv_function, groups), total=len(groups),
                         desc="Calculating n-values: ",
                         leave=True))

        elif settings.debug_level >= 1:
            for group in tqdm(groups, desc="Calculating n-values: ", total=len(groups)):
                data, n = NValueCalculator.analyze_group(self, group, settings.save_n_value_data, settings.graph_n_value_calculations)
                results.append(data)
                n_data.append(n)

            for x in n_data:
                if x:
                    n_value_data.append(x[0])

        if settings.save_n_value_data:
            prepared_df = pd.concat(results)
            n_value_df = pd.DataFrame(n_value_data, columns=["Index", "cf", "sample_group", "n_value", "std_dev",
                                                             "all_std_dev", "all_n_D_values"])
            n_value_df.to_csv(self.out_path[:-35] + "n_value_data.tsv", sep='\t')
        else:
            prepared_df = pd.concat(results)

        prepared_df = prepared_df.set_index('index')
        prepared_df = prepared_df.sort_index()

        self.full_df = self.full_df.merge(right=prepared_df[output_columns],
                                          how='outer',
                                          left_index=True,
                                          right_index=True)

    def analyze_group(self, partition, save_data, make_graphs):
        """
        Generates emass information for each chemical formula and then sends each possible charge state to analyze_row()

        Parameters:
            partition (pd.DataFrameGroupBy): A DataFrameGroupBy that contains data for a single chemical formula. Can contain
                multiple rows, each representing a different isotopic envelope.
            save_data (bool): whether we will create n-value calculation data .tsv output
            make_graphs (bool): whether we will create graphs showing n-value calculation process

        Returns:
            pd.DataFrame: A DataFrame containing the n-value and associated standard deviation and detailed n-value calculation
            data if the appropriate setting is turned on.
        """

        results = []
        nvalues = []
        filtered_out_nvalues = []
        reason_for_filtering = []

        try:
            # Each chemical formula has only 1 enrichment value Because adducts can occur,
            # we need to only look at hydrogen that would be in the actual molecule. I need to ask Rusty if I should
            # just choose the smaller cf, or which cf we want to use.
            if self.biomolecule_type == "Peptide":
                num_h_adduct = np.inf
                num_h_no_adduct, chem_f, _ = NValueCalculator.parse_cf(partition[1].iloc[0]["chemical_formula"])
            else:
                num_h_adduct, chem_f, _ = NValueCalculator.parse_cf(partition[1].iloc[0]["adduct_chemical_formula"])
                num_h_no_adduct, _, _ = NValueCalculator.parse_cf(partition[1].iloc[0]["chemical_formula"])

            if num_h_no_adduct < num_h_adduct:
                num_h = num_h_no_adduct
            else:
                num_h = num_h_adduct

            n_begin = 0
            low_pct = 0
            num_peaks = max(partition[1].intensities.str.len())
            emass_results_dict = dict()

            # sort the rows, so we have the rows to calculate n-values for at the top
            partition[1].sort_values(by='calculate_n_value', ascending=False, kind='mergesort', ignore_index=True,
                                     inplace=True)

            # calculate n-values for appropriate rows, once we get to a row that has 'no' in the calculate_n_value column,
            # we break out of the for loop and average the n-values we have
            n_value_data = []
            for row in partition[1].itertuples(index=False):
                if row.calculate_n_value.lower() == 'no':
                    break
                # Get/find enrichment value
                enrichment = row.enrichment
                if enrichment == 0:
                    enrichment = 0.05
                if enrichment not in emass_results_dict:
                    high_pct = enrichment
                    try:
                        # set n_begin to zero, Coleman 2025. This is a distribution.
                        emass_results = emass(chem_f, 0, num_h, low_pct, high_pct, num_peaks)
                        emass_results_dict[enrichment] = emass_results
                    except ValueError as e:
                        print("Chemical Formula ", row.cf, " contains unsupported molecules")
                emass_results = emass_results_dict[enrichment]

                # Calculate n-value and standard deviation, gather n-value calculation data if setting is turned on
                try:
                    n_value, stddev, n_data, med_first_derivative, time_point = NValueCalculator.analyze_row(self, row,
                                                                                                             emass_results, save_data, make_graphs)

                    # Only allow the append if the stddev is less than .02
                    # if stddev <= 0.02: #make this a setting? Coleman 2025
                    if med_first_derivative >= 0.2 and stddev <= 0.05:

                        nvalues.append([n_value, stddev, time_point]) # remove time_point from here if it does not work - Coleman
                    else:
                        # Determine reason(s) for filtering
                        reasons = []
                        if med_first_derivative < 0.2:
                            reasons.append("derivative too small")
                        if stddev > 0.05:
                            reasons.append("stddev too large")

                        # Append to filtered list and reason list
                        filtered_out_nvalues.append([n_value, stddev, time_point])
                        reason_for_filtering.append(" and ".join(reasons))

                    if save_data:
                        n_value_data.append(n_data)

                except ValueError as e:
                    print("EXCEPTION OCCURRED WITH {}!".format(row))

            has_error = False
            nv_data = None

            # append error message if no valid time points could be used
            if len(nvalues) == 0:
                has_error = True
                data = timepoints_statement

            if not has_error:
                # calculate standard deviation of n values and exclude any outliers - Ben D
                nv = [n[0] if n[0] is not np.nan else -1 for n in nvalues]
                FN_stddev = [n[1] if n[0] is not np.nan else -1 for n in nvalues]
                time_points = [n[2] if n[0] is not np.nan else -1 for n in nvalues]
                filtered_out_nv = [n[0] if n[0] is not np.nan else -1 for n in filtered_out_nvalues]
                filtered_out_FN_stddev = [n[1] if n[0] is not np.nan else -1 for n in filtered_out_nvalues]
                filtered_out_times = [n[2] if n[0] is not np.nan else -1 for n in filtered_out_nvalues]

                nv_std = np.std(nv, dtype=np.float64)
                avg_nv = np.median(nv)  # change this to median, Coleman
                avg_FN_std = np.median(FN_stddev)
                mad_nv = np.median(np.abs(nv - avg_nv))  # MAD about the median

                # Guard against the (rare) case where all nv are identical
                if mad_nv == 0:
                    mad_nv = np.finfo(float).eps  # smallest positive float

                # below, instead of removing outliers, we use the median
                '''
                for i, value in enumerate(nv):
                    # remove any n values that aren't within 1 standard deviation - Ben D   Return this if necessary, Coleman 
                    if value > avg_nv + (2*nv_std) or value < avg_nv - (2* nv_std):
                        filtered_out_nv.append(nv.pop(i))
                        filtered_out_FN_stddev.append(FN_stddev.pop(i))
                        filtered_out_times.append(time_points.pop(i))
                        reason_for_filtering.append('outside 2 stddev')
                '''

                threshold = 3.5
                scale = 0.6745 * mad_nv  # rescales MAD to σ‑equivalent

                # iterate *backwards* so list pops don’t shift later indices
                for i in range(len(nv) - 1, -1, -1):
                    modified_z = (nv[i] - avg_nv) / scale
                    if abs(modified_z) > threshold:
                        filtered_out_nv.append(nv.pop(i))
                        filtered_out_FN_stddev.append(FN_stddev.pop(i))
                        filtered_out_times.append(time_points.pop(i))
                        reason_for_filtering.append('|modified z| > 3.5 (MAD)')

                # if there aren't any values left, append error message
                if not nv:
                    has_error = True
                    nv_data = timepoints_statement

                if not has_error:
                    # recalculate average n-value and standard deviation with outliers removed - Ben Driggs
                    avg_nv = np.median(nv)
                    nv_std = np.std(nv, dtype=np.float64)
                    num_points = len(nv)
                    avg_FN_std = np.median(FN_stddev)
                    lower_margin, upper_margin = bootstrap_median_ci(nv)  # remove this if it breaks things coleman 2025

                    # calculate cv (confidence value - stddev/n-value)
                    cv = float(nv_std / avg_nv)

                    # added the median of the FN_std
                    nv_data = [float(avg_nv), float(nv_std), int(num_points), float(cv), float(avg_FN_std), float(lower_margin),
                               float(upper_margin), nv, time_points, filtered_out_nv, filtered_out_times, reason_for_filtering]

            # apply average n value to each time point - Ben D
            for row in partition[1].itertuples(index=False):
                try:
                    results.append(nv_data)
                except Exception as e:
                    print(e)
                    print("EXCEPTION OCCURRED WITH {}!".format(row))

            if save_data or make_graphs:
                return pd.concat([partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
                                 axis=1), n_value_data
            else:
                return pd.concat([partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
                                 axis=1)
        except IOError as e:
            print(e)
            return pd.concat(
                [partition[1].reset_index(drop=True), pd.DataFrame(data=error_statement, columns=output_columns)],
                axis=1), []

    def analyze_row(self, row, emass_results, save_data, make_graphs):
        """
        Calculates n-value for each row

        Parameters:
            row (pd.Series): Contains chemical formula, empirical intensity data, and enrichment level, and whether an n-value should be calculated
            emass_results (pd.DataFrame): Results from emass containing unlabeled & labeled intensity and m/z data.
            save_data (bool): whether we will create n-value calculation data .tsv output
            make_graphs (bool): whether we will create graphs showing n-value calculation process

        Returns:
            Tuple: a tuple containing the given n-value, the stddev associated with it (will be empty if n-value was not calculated), and
            detailed info about the n-value calculation process if the appropriate setting is turned on.
        """

        def dIt_filter(unfiltered_dataframe, dIt_data):
            return_value = unfiltered_dataframe.where(dIt_data > 0.05, np.nan)
            return_value['n_D'][0] = 0.0
            return return_value

        def n_value_by_stddev(fraction_new, MINIMUM_PEAKS=3):
            try:
                # Compute standard deviation dataframe
                stddev_dataframe = calc_stddev(fraction_new, MINIMUM_PEAKS)

                # Ensure n_D values are finite and exclude NaNs
                valid_entries = stddev_dataframe  # temp Coleman
                # valid_entries = stddev_dataframe[['n_D', 'n_value_stddev']].dropna()
                # valid_entries = valid_entries[np.isfinite(valid_entries['n_D'])]

                if valid_entries.empty:
                    filtered_nValue_min = np.nan
                    filtered_stddev_min = np.nan

                    truth_values = ~np.isnan(stddev_dataframe.iloc[:, 2:-1])
                    truth_count = [np.sum(truth_values.iloc[:, x]) for x in range(truth_values.shape[1])]
                    peaks_included = ' '.join(str(x) for x in truth_count)
                else:
                    min_index = valid_entries['n_value_corrected_stddev'].idxmin()

                    filtered_nValue_min = valid_entries.iloc[min_index]['n_D']
                    filtered_stddev_min = valid_entries.iloc[min_index]['n_value_corrected_stddev']

                    data_mask = stddev_dataframe.iloc[min_index].notna().reset_index(drop=True)[2:-1]
                    peaks_included = pd.Series(stddev_dataframe.columns)[2:-1][data_mask] \
                        .to_string(header=False, index=False).replace('I', 'M').replace('\n', ' ')
                    if peaks_included == 'Series([], )':
                        peaks_included = ''

                return stddev_dataframe, filtered_nValue_min, filtered_stddev_min, peaks_included

            except Exception as e:
                print(f"An error occurred on line {e.__traceback__.tb_lineno}: {e}")
                traceback.print_exc()
                raise

        def n_value_dIt_filter(unfiltered_fraction_new, dIt_data, MINIMUM_PEAKS=3):
            # Generate dIt filtered fraction_new DataFrame
            filtered_fraction_new = dIt_filter(unfiltered_fraction_new, dIt_data)

            try:
                filtered_fraction_new.drop('n_value_stddev', axis=1, inplace=True)
                filtered_fraction_new.drop('num_non_null', axis=1, inplace=True)
            except:
                pass
            return n_value_by_stddev(filtered_fraction_new)

        def calc_stddev(fraction_new, MINIMUM_PEAKS=3):
            stddev_dataframe = fraction_new.copy()

            # Find how many valid peaks exist
            stddev_dataframe = stddev_dataframe[['n_D'] + [col for col in stddev_dataframe.columns if 'I' in col]]
            stddev_dataframe['num_non_null'] = (stddev_dataframe.count(axis='columns') - 1)

            # rearrange columns to allow for proper std_dev calculation
            stddev_dataframe = stddev_dataframe[
                ['num_non_null'] + [col for col in stddev_dataframe.columns if col != 'num_non_null']]

            # Calculate stddev for valid rows
            stddev_dataframe['n_value_stddev'] = stddev_dataframe.loc[:, 'I0'::].std(axis='columns')
            stddev_dataframe['n_value_corrected_stddev'] = stddev_dataframe['n_value_stddev'] * stddev_dataframe['n_D']

            # used to be >= MINIMUM_PEAKS, but changed it to >= MINIMUM_PEAKS-1 to allow calculating std dev for rows with only 2 peaks
            # consider changing this if it hurts results/statistics - Ben D
            stddev_dataframe.loc[:, 'n_value_stddev'] = stddev_dataframe['n_value_stddev'].where(
                stddev_dataframe['num_non_null'] >= MINIMUM_PEAKS, np.nan)

            return stddev_dataframe

        def s_n_filter(empirical_intensities, noise, NOISE_FILTER=10.0):
            empir_series = pd.Series(empirical_intensities)
            s_n_values = pd.Series([empirical_intensities[i] / noise for i in range(len(empirical_intensities))])
            filtered_empirical_intensities = empir_series.where(s_n_values >= NOISE_FILTER, np.nan)
            return filtered_empirical_intensities

        def dIe_filter(unfiltered_dataframe, dIe_data):
            filtered_dataframe = unfiltered_dataframe.where(dIe_data < 0.05, np.nan)
            filtered_dataframe['n_D'] = unfiltered_dataframe['n_D']
            return filtered_dataframe

        def n_value_dIe_filter(unfiltered_fraction_new, dIe_data, MINIMUM_PEAKS=3):
            # Generate dIt filtered fraction_new DataFrame
            filtered_fraction_new = dIe_filter(unfiltered_fraction_new, dIe_data)
            filtered_fraction_new['n_D'] = unfiltered_fraction_new['n_D']
            filtered_fraction_new.drop('n_value_stddev', axis=1, inplace=True)
            filtered_fraction_new.drop('num_non_null', axis=1, inplace=True)

            return n_value_by_stddev(filtered_fraction_new)

        def n_value_s_n_filter(empirical_intensities, noise, fraction_new_values, NOISE_FILTER=100.0):
            valid_intensities = s_n_filter(empirical_intensities, noise, NOISE_FILTER)
            filtered_fraction_new = fraction_new_values.copy()
            for p, i in enumerate(valid_intensities):
                if np.isnan(i):
                    filtered_fraction_new['I' + str(p)] = np.nan
            return n_value_by_stddev(filtered_fraction_new)

        def n_value_using_angles(empirical_intensities, emass_labeled_intensities, MIN_DIMENSION=2):
            def angle_between(v1, v2):
                """Calculates the unit vector of a given vector
                Parameters
                ----------
                v1 : :obj:`list` of :obj:`float`
                    The first vector to compare
                v2 : :obj:`list` of :obj:`float`
                    The second vector to compare
                Returns
                ----------
                :obj:`float`
                    Angle between the two vectors, given in radians
                """

                def unit_vector(vector):
                    """Calculates the unit vector of a given vector
                    Parameters
                    ----------
                    vector : :obj:`list` of :obj:`float`
                    Returns
                    ----------
                    vector : :obj:`list` of :obj:`float`
                    """
                    return vector / np.linalg.norm(vector)

                if len(v1) <= MIN_DIMENSION:
                    return np.nan

                uv1 = unit_vector(v1)
                uv2 = unit_vector(v2)
                angle = np.arccos(np.dot(uv1, uv2))
                if np.isnan(angle):
                    if (uv1 == uv2).all():
                        return 0.0
                    else:
                        return np.pi
                return angle

            # Turn empirical intensities into a vector to compare angles.
            empirical_vector = np.array(empirical_intensities)

            # Turn emass outputs into vectors to use to compare angles
            emass_labeled_vectors = emass_labeled_intensities.drop(columns="n_D")
            emass_labeled_vectors['combined_intensities'] = emass_labeled_vectors.values.tolist()
            emass_labeled_vectors['combined_intensities'] = emass_labeled_vectors['combined_intensities'].apply(
                lambda x: np.array(x))
            emass_labeled_vectors = emass_labeled_vectors['combined_intensities']

            empir_truths = ~np.isnan(empirical_vector)
            total_truths = np.array([~np.isnan(emass_labeled_vectors[i]) & empir_truths for i in
                                     np.array(range(len(emass_labeled_vectors)))])

            # Compute angle between emass intensity data and empirical intensity data
            empir_v_emass_angles = pd.Series(
                [angle_between(empirical_vector[total_truths[i]], emass_labeled_vectors[i][total_truths[i]]) for i in
                 range(len(emass_labeled_vectors))])

            # Determine n-value (the #D associated with the smallest angle)
            angle_n_value = empir_v_emass_angles.to_numpy().argmin()

            column_names = emass_labeled_intensities.columns[1:]
            peaks_included = pd.Series(column_names[total_truths[angle_n_value]]) \
                .to_string(header=False, index=False).replace('I', 'M').replace('\n', ' ')
            if peaks_included == 'Series([], )':
                truth_count = [np.sum(total_truths[:, x]) for x in range(len(total_truths[0]))]
                peaks_included = ' '.join(str(x) for x in truth_count)

            return empir_v_emass_angles, angle_n_value, peaks_included

        # could speed up by forgoing mzs, but mzs will be needed down the line
        emass_unlabeled_data, emass_labeled_data = emass_results

        def interpolate_mz_values(df):
            interpolated_rows = []

            try:
                # Ensure the original dataframe's n_D column is converted to float (handles string integers)
                df['n_D'] = df['n_D'].astype(float)
            except Exception as e:
                print(f"Error converting 'n_D' to float: {e}")
                return df  # Return original DataFrame if conversion fails

            try:
                # Loop through each pair of adjacent n_D values
                for i in range(len(df) - 1):
                    try:
                        # Accessing by column names (strings), not by index positions
                        n1, n2 = df.iloc[i]['n_D'], df.iloc[i + 1]['n_D']
                        mz_0_1 = df.iloc[i]['mz0']  # mz_0 stays the same
                        mz_1_1, mz_1_2 = df.iloc[i]['mz1'], df.iloc[i + 1]['mz1']
                        mz_2_1, mz_2_2 = df.iloc[i]['mz2'], df.iloc[i + 1]['mz2']
                    except KeyError as e:
                        print(f"Missing column: {e}")
                        return df
                    except Exception as e:
                        print(f"Error accessing row {i}: {e}")
                        return df

                    # Create interpolated values for n = n1.1 to n = n2.9
                    for j in range(1, 10):  # for n = n1.1 to n1.9, then n2.0 to n2.9
                        try:
                            alpha = j / 10.0  # alpha goes from 0.1 to 0.9
                            n_interpolated = n1 + alpha

                            # Keep mz_0 constant and interpolate mz_1 and mz_2
                            mz_0_interpolated = mz_0_1
                            mz_1_interpolated = (1 - alpha) * mz_1_1 + alpha * mz_1_2
                            mz_2_interpolated = (1 - alpha) * mz_2_1 + alpha * mz_2_2

                            interpolated_rows.append(
                                [n_interpolated, mz_0_interpolated, mz_1_interpolated, mz_2_interpolated])
                        except Exception as e:
                            print(f"Error interpolating at n_D={n1} (step {j}): {e}")
                            continue  # Skip this interpolation step if it fails

            except Exception as e:
                print(f"Unexpected error in interpolation loop: {e}")
                return df

            try:
                # Convert the list of interpolated rows into a new DataFrame
                interpolated_df = pd.DataFrame(interpolated_rows, columns=['n_D', 'mz0', 'mz1', 'mz2'])

                # Concatenate and sort
                result_df = pd.concat([df, interpolated_df], ignore_index=True).sort_values(by='n_D').reset_index(
                    drop=True)
                return result_df
            except Exception as e:
                print(f"Error creating the final DataFrame: {e}")
                return df

        def interpolate_n_values(df):
            interpolated_rows = []

            try:
                # Loop through each pair of adjacent n_D values
                for i in range(len(df) - 1):
                    try:
                        # Accessing by column names (strings), not by index positions
                        n1, n2 = df.iloc[i]['n_D'], df.iloc[i + 1]['n_D']
                        I_0_1, I_0_2 = df.iloc[i]['I0'], df.iloc[i + 1]['I0']
                        I_1_1, I_1_2 = df.iloc[i]['I1'], df.iloc[i + 1]['I1']
                        I_2_1, I_2_2 = df.iloc[i]['I2'], df.iloc[i + 1]['I2']
                    except KeyError as e:
                        print(f"Missing column: {e}")
                        return df
                    except Exception as e:
                        print(f"Error accessing row {i}: {e}")
                        return df

                    # Create interpolated values for n = n1.1 to n = n2.9
                    for j in range(1, 10):  # for n = n1.1 to n1.9, then n2.0 to n2.9
                        try:
                            alpha = j / 10.0  # alpha goes from 0.1 to 0.9
                            n_interpolated = n1 + alpha

                            I_0_interpolated = (1 - alpha) * I_0_1 + alpha * I_0_2
                            I_1_interpolated = (1 - alpha) * I_1_1 + alpha * I_1_2
                            I_2_interpolated = (1 - alpha) * I_2_1 + alpha * I_2_2

                            interpolated_rows.append(
                                [n_interpolated, I_0_interpolated, I_1_interpolated, I_2_interpolated])
                        except Exception as e:
                            print(f"Error interpolating at n_D={n1} (step {j}): {e}")
                            continue  # Skip this interpolation step if it fails

            except Exception as e:
                print(f"Unexpected error in interpolation loop: {e}")
                return df

            try:
                # Convert the list of interpolated rows into a new DataFrame
                interpolated_df = pd.DataFrame(interpolated_rows, columns=['n_D', 'I0', 'I1', 'I2'])

                # Concatenate and sort
                result_df = pd.concat([df, interpolated_df], ignore_index=True).sort_values(by='n_D').reset_index(
                    drop=True)
                return result_df
            except Exception as e:
                print(f"Error creating the final DataFrame: {e}")
                return df

        # Returns 2 tuples of DataFrame with rows representing # of hydrogen that have been deuterated,
        # columns represent the peak #, data is m/z for 1st and intensity for 2nd
        (emass_unlabeled_mz, emass_unlabeled_intensities) = emass_unlabeled_data
        (emass_labeled_mz, emass_labeled_intensities) = emass_labeled_data

        if settings.interpolate_n_values:
            emass_labeled_intensities = interpolate_n_values(emass_labeled_intensities)

        # Get the length of emass_labeled_data
        target_length = len(emass_labeled_intensities)

        # Repeat the first row of emass_unlabeled_data to match the target length
        if not emass_unlabeled_intensities.empty:
            first_row = emass_unlabeled_intensities.iloc[0]  # Extract first row
            emass_unlabeled_intensities = pd.DataFrame([first_row] * target_length).reset_index(drop=True)
            emass_unlabeled_intensities["n_D"] = emass_labeled_intensities['n_D']

        # Repeat for emass_unlabeled_mz if needed
        if not emass_unlabeled_mz.empty:
            first_row = emass_unlabeled_mz.iloc[0]  # Extract first row
            emass_unlabeled_mz = pd.DataFrame([first_row] * target_length).reset_index(drop=True)
            emass_unlabeled_mz["n_D"] = emass_labeled_intensities['n_D']

        unfiltered_fraction_new = emass_unlabeled_intensities.copy()
        dIt_data = emass_unlabeled_intensities.copy()

        # Generate fraction_new and dIt data
        for peak, ie in enumerate(row.intensities):
            # calculate fraction new for theoretical emass distributions
            unfiltered_fraction_new['I' + str(peak)] = (
                    (emass_unlabeled_intensities['I' + str(peak)].iloc[0] - ie) /
                    (emass_unlabeled_intensities['I' + str(peak)].iloc[0] - emass_labeled_intensities['I' + str(peak)])
            )

            # calculates deltas for theoretical emass distributions, later used in n_value_dIt_filter to remove peaks
            # that have deltas less than 0.05
            dIt_data['I' + str(peak)] = (
                np.abs((emass_labeled_intensities['I' + str(peak)] - emass_unlabeled_intensities['I' + str(peak)].iloc[
                    0])))

        # Calculate n-value with no filters (using stddev)
        unfiltered_fraction_new, n_value, stddev, peaks_included = n_value_by_stddev(unfiltered_fraction_new)

        # Define the domain limit
        n_min = n_value - 10
        n_max = n_value + 10

        # Filter the data to the desired range
        filtered_data = unfiltered_fraction_new[
            (unfiltered_fraction_new['n_D'] >= n_min) &
            (unfiltered_fraction_new['n_D'] <= n_max)
            ]

        # Calculate first derivative and its absolute value
        first_derivative = abs(np.gradient(filtered_data['n_value_corrected_stddev'], filtered_data['n_D']))
        med_first_derivative = np.median(first_derivative)  # Coleman 2025

        # filters out peaks that have a delta of less than 0.05
        dIt_unfiltered_fraction_new, dIt_n_value, dIt_stddev, dIt_peaks_included = n_value_dIt_filter(unfiltered_fraction_new, dIt_data)

        def sanitize_filename(string):
            """Removes invalid characters from filenames."""
            return re.sub(r'[<>:"/\\|?*]', '-', string)

            # if setting is turned on, we'll graph out a visualization of how we calculate n-values
        n_value_info = []

        # if the make_graphs and stddev < 0.02: #remove the 0.2 filter for now
        if med_first_derivative >= 0.2 and stddev <= 0.05:  # are we being too conservative here? 3/12/2025 -coleman
            if not np.isnan(n_value):
                # Plotting the unfiltered data
                plt.plot(unfiltered_fraction_new['n_D'], unfiltered_fraction_new['n_value_corrected_stddev'],
                         marker='.', linestyle='none', label="dIt n")

                # Plot first derivative
                plt.plot(filtered_data['n_D'], first_derivative, label="First Derivative", color='orange')

                # Plot a horizontal line at the average first derivative value
                plt.axhline(y=med_first_derivative, color='blue', linestyle='--', label="Med First Derivative")

                # Title and labels
                if self.biomolecule_type == 'Lipid':
                    plt.title(str.format("{}_{}\n n_value = {}, std_dev = {},\n Median First Derivative = {:.2f}",
                                         row._1, row.bio_rep, n_value, stddev, med_first_derivative))
                else:
                    # changed avg_first_derivative to med_first_derivative, if this breaks this is probably why, Coleman
                    plt.title(str.format("{}_{}\n dIt_n_value = {}, dIt_std_dev = {}, Median First Derivative = {:.2f}",
                                         row.chemical_formula, row.sample_group, int(dIt_n_value), int(dIt_stddev),
                                         med_first_derivative))

                plt.xlabel("Number of Deuteriums")
                plt.ylabel("Standard Deviation (corrected for n-value)")
                plt.xlim(n_value - 30, n_value + 30)
                plt.ylim(0, 10)
                plt.plot(n_value, stddev, marker="o", markeredgecolor="red")
                plt.legend()

                # Save the figure
                name = sanitize_filename(row._1)
                plt.savefig(str.format(self.graph_folder + "\\{}_{}.png", name, row.bio_rep))
                plt.clf()

        # send back n-value calculation info if the setting is turned on, otherwise we'll just return an empty list and move on
        if save_data:
            n_value_info = [row.index, row.chemical_formula, row.sample_group, n_value,
                            stddev, unfiltered_fraction_new['n_value_corrected_stddev'].values,
                            unfiltered_fraction_new['n_D'].values]

        return n_value, stddev, n_value_info, med_first_derivative, row.time

    @staticmethod
    def parse_cf(chem_f):
        d = dict(re.findall(r'([A-Z][a-z]*)(\d*)', chem_f))

        num_h = int(d.get('H', '0'))  # Default to 0 if 'H' is not present
        num_d = int(d.get('D', '0'))  # Default to 0 if 'D' is not present

        # Add Deuterium count to Hydrogen and remove 'D' entry
        num_h += num_d
        d.pop('D', None)  # Remove 'D' safely if it exists

        # Replace Hydrogen count with '{}' and add placeholder for labeled sites
        d['H'] = '{}'
        d['X'] = '{}'

        cf_string = ''.join(f"{k}{v}" for k, v in d.items())

        return num_h, cf_string, num_d


def main():
    filename = "../n_value_debug_data.tsv"
    settings_path = "../resources/temp_settings.yaml"
    df = pd.read_csv(filename, sep='\t')
    calculator = NValueCalculator(df, settings_path, "Lipid", "../n_value_test")
    calculator.run()
    calculator.full_df.to_csv(filename[:-4] + "_nvalue.tsv", sep='\t')


if __name__ == "__main__":
    main()
