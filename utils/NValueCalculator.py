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
except:
    from n_value_calc_emass import emass

from functools import partial
import pandas as pd
import re
import time
import multiprocessing as mp
import concurrent.futures as cf
import numpy as np
from tqdm import tqdm
import deuterater.settings as settings

# output_columns = ['empir_n', 'stddev', 'dIt_n', 'dIt_stddev']
output_columns = ['n_value', 'n_value_stddev', "num_nv_time_points", "cv"]
error_statement = ["error occurred", "error occurred", "error occurred", "error occurred"]
timepoints_statement = ["no valid time points", "no valid time points", "no valid time points", "no valid time points"]


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

    def prepare_model(self):
        if "Adduct_cf" in self.full_df.columns:
            df = self.full_df.rename(columns={  # TODO: Verify Column Names
                'Adduct_cf': 'adduct_chemical_formula',
                'cf': 'chemical_formula',
                'abundances': 'intensities',
                'labeling': 'enrichment'})
        else:
            # Renames columns to proper names depending on if a lipid or protein is being passed in
            df = self.full_df.rename(columns={  # TODO: Verify Column Names
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
        df = df[['index', 'chemical_formula', 'adduct_chemical_formula', 'intensities', 'enrichment', 'sample_group',
                 'bio_rep', 'calculate_n_value']]

        # TODO: should we divide by genotype (sample_group), not bio_rep? - Ben Driggs
        if self.biomolecule_type == "Peptide":
            df.loc[:, "divider"] = df['chemical_formula'] + df['sample_group']
        else:
            df.loc[:, "divider"] = df['chemical_formula'] + df['adduct_chemical_formula'] + df['sample_group']
        df.sort_values(by="divider")

        self.prepared_df = df

    def run(self):
        self.prepare_model()

        groups = self.prepared_df.groupby('divider', sort=False)

        nv_function = partial(NValueCalculator.analyze_group, self, save_data=settings.save_n_value_data, make_graphs=settings.graph_n_value_calculations)

        results = list()
        n_value_data = list()
        
        # if either of these settings are turned on, we need to run this step through debug mode so we can get all the results back
        if settings.graph_n_value_calculations or settings.save_n_value_data:
            settings.debug_level = 1

        settings.debug_level = 1
        n_data = []
        if settings.debug_level == 0:
            with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                results = list(
                    tqdm(executor.map(nv_function, groups), total=len(groups),
                         desc="Calculating n-values: ",
                         leave=True))

        elif settings.debug_level >= 1:
            for group in tqdm(groups, desc="Calculating n-values: ", total=len(groups)):
                data = NValueCalculator.analyze_group(self, group, settings.save_n_value_data, settings.graph_n_value_calculations)
                results.append(data[0])
                n_data.append(data[1])
                
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

        try:
            # print(f"Started group #" + str(partition[1].index[0]) + f" at: " + str(time.perf_counter() - start) +
            # f" seconds", end='\r') Each chemical formula has only 1 enrichment value Because adducts can occur,
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

            # num_h, cf, _ = NValueCalculator.parse_cf(partition[0])
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
                        emass_results = emass(chem_f, n_begin, num_h, low_pct, high_pct, num_peaks)
                        emass_results_dict[enrichment] = emass_results
                    except Exception as e:
                        try:
                            print("Chemical Formula ", row.cf, " contains unsupported molecules")
                        finally:
                            print("Chemical Formula ", row.chemical_formula, " contains unsupported molecules")
                emass_results = emass_results_dict[enrichment]

                # Calculate n-value and standard deviation, gather n-value calculation data if setting is turned on
                try:
                    n_value, stddev, n_data = NValueCalculator.analyze_row(self, row, emass_results, save_data, make_graphs)
                    nvalues.append([n_value, stddev])

                    if save_data:
                        n_value_data.append(n_data)

                except Exception as e:
                    print("EXCEPTION OCCURRED WITH {}!".format(row))

            has_error = False

            # append error message if no valid time points could be used
            if len(nvalues) == 0:
                has_error = True
                data = timepoints_statement

            if not has_error:
                # calculate standard deviation of n values and exclude any outliers - Ben D
                nv = [n[0] if n[0] is not np.nan else -1 for n in nvalues]
                nv_std = np.std(nv, dtype=np.float64)
                avg_nv = np.average(nv)

                for i, value in enumerate(nv):
                    # remove any n values that aren't within 1 standard deviation - Ben D
                    if value > avg_nv + nv_std or value < avg_nv - nv_std:
                        nv.pop(i)

                # if there aren't any values left, append error message
                if not nv:
                    has_error = True
                    data = timepoints_statement

                if not has_error:
                    # recalculate average n-value and standard deviation with outliers removed - Ben Driggs
                    avg_nv = np.average(nv)
                    nv_std = np.std(nv, dtype=np.float64)
                    num_points = len(nv)

                    # calculate cv (confidence value - stddev/n-value)
                    cv = float(nv_std / avg_nv)
                    data = [float(avg_nv), float(nv_std), int(num_points), float(cv)]

            # apply average n value to each time point - Ben D
            for row in partition[1].itertuples(index=False):
                try:
                    results.append(data)
                except Exception as e:
                    print(e)
                    print("EXCEPTION OCCURRED WITH {}!".format(row))

            if save_data or make_graphs:
                return pd.concat([partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
                                 axis=1), n_value_data
            else:
                return pd.concat([partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
                                 axis=1), []
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
            stddev_dataframe = calc_stddev(fraction_new, MINIMUM_PEAKS)

            # Find n-value and stddev of filtered data
            if stddev_dataframe['n_value_stddev'].isna().all():
                filtered_nValue_min = np.nan
                filtered_stddev_min = np.nan

                truth_values = ~np.isnan(stddev_dataframe.iloc[:, 2:-1])
                truth_count = [np.sum(truth_values.iloc[:, x]) for x in range(truth_values.shape[1])]
                peaks_included = ' '.join(str(x) for x in truth_count)
            else:
                min_index = stddev_dataframe['n_value_stddev'].idxmin()
                filtered_nValue_min = stddev_dataframe.loc[min_index, 'n_D']
                filtered_stddev_min = stddev_dataframe['n_value_stddev'].min()

                data_mask = stddev_dataframe.iloc[min_index].notna().reset_index(drop=True)[2:-1]
                peaks_included = pd.Series(stddev_dataframe.columns)[2:-1][data_mask] \
                    .to_string(header=False, index=False).replace('I', 'M').replace('\n', ' ')
                if peaks_included == 'Series([], )':
                    peaks_included = ''
            return stddev_dataframe, filtered_nValue_min, filtered_stddev_min, peaks_included

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
            for peak, ie in enumerate(valid_intensities):
                if (np.isnan(ie)):
                    filtered_fraction_new['I' + str(peak)] = np.nan
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

        # Returns 2 tuples of DataFrame with rows representing # of hydrogen that have been deuterated,
        # columns represent the peak #, data is m/z for 1st and intensity for 2nd
        (emass_unlabeled_mz, emass_unlabeled_intensities) = emass_unlabeled_data
        (emass_labeled_mz, emass_labeled_intensities) = emass_labeled_data

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

        # filters out peaks that have a delta of less than 0.05
        dIt_unfiltered_fraction_new, dIt_n_value, dIt_stddev, dIt_peaks_included = n_value_dIt_filter(
            unfiltered_fraction_new, dIt_data)

        # if setting is turned on, we'll graph out a visualization of how we calculate n-values
        n_value_info = []
        if make_graphs:
            import matplotlib.pyplot as plt
            if not np.isnan(dIt_n_value):
                plt.plot(dIt_unfiltered_fraction_new['n_value_stddev'], label="dIt n")
                plt.title(
                    str.format("{}_{}\n dIt_n_value = {}, dIt_std_dev = {}", row.chemical_formula, row.sample_group,
                               int(dIt_n_value), int(dIt_stddev)))
                plt.legend()
                plt.xlabel("Number of Deuteriums")
                plt.ylabel("Standard Deviation")
                plt.plot(dIt_n_value, dIt_stddev, marker="o", markeredgecolor="red")
                plt.savefig(str.format(self.graph_folder + "\\{}_{}_dIt_{}.png", row.chemical_formula, row.sample_group,
                                       row.index))
                plt.clf()
                
        # send back n-value calculation info if the setting is turned on, otherwise we'll just return an empty list and move on
        if save_data:
            n_value_info = [row.index, row.chemical_formula, row.sample_group, dIt_n_value,
                            dIt_stddev, dIt_unfiltered_fraction_new['n_value_stddev'].values,
                            dIt_unfiltered_fraction_new['n_D'].values]

        return dIt_n_value, dIt_stddev, n_value_info

    @staticmethod
    def parse_cf(cf):
        d = dict(re.findall(r'([A-Z][a-z]*)(\d*)', cf))
        num_h = int(d['H'])

        num_d = 0
        # Makes the cf not distinguish between Deuterium and Hydrogen
        # if d['D']:
        # 	num_h += int(d['D'])
        # 	del d['D']

        blank = '{}'
        d.update({'H': blank, 'X': blank})  # H = # of Hydrogen, X = # labeled sites
        cf_string = ''.join('%s%s' % (k, v) for k, v in d.items())
        return num_h, cf_string, num_d


def main():
    filename = "../n_value_debug_data.tsv"
    settings_path = "../resources/temp_settings.yaml"
    df = pd.read_csv(filename, sep='\t')
    calculator = NValueCalculator(df, settings_path, "Lipid")
    calculator.run()
    calculator.full_df.to_csv(filename[:-4] + "_nvalue.tsv", sep='\t')


if __name__ == "__main__":
    main()
