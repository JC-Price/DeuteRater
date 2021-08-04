from utils.n_value_calc_emass import emass
import pandas as pd
import re
import time
import multiprocessing as mp
import numpy as np

start = time.perf_counter()

output_columns = ['empir_n', 'stddev']

class NValueCalculator:
	def __init__(self, dataframe, output_columns=['']):
		# this DataFrame should only contain the correct columns
		self.df = dataframe
		self.output_columns = output_columns

	def run(self):
		groups = self.df.groupby('chemical_formula', sort=False)

		try:
			cpus = mp.cpu_count()
		except:
			cpus = 2  # Default value

		# results = list()
		# for group in groups:
		# 	results.append(self.analyze_group(group))

		with mp.Pool(processes=cpus) as pool:
			from itertools import product
			toIter = product(groups)
			results = pool.starmap(NValueCalculator.analyze_group, toIter)

		self.df = pd.concat(results)

	@staticmethod
	def analyze_group(partition):
		"""
        Generates emass information for each chemical formula and then sends each possible charge state to analyze_row()

        Parameters:
            partition (pd.DataFrameGroupBy): A DataFrameGroupBy that contains data for a single chemical formula. Can contain
                multiple rows, each representing a different isotopic envelope.
            plot_dir (str): The directory to output fraction_new .csv files to.
            enrichment (float): The %D enrichment level of the given partition.

        Returns:
            pd.DataFrame: A DataFrame containing the n-value and associated standard deviation
        """
		results = []
		try:
			print(f"Started group #" + str(partition[1].index[0]) + f" at: " + str(time.perf_counter() - start) + f" seconds", end='\r')
			# Each chemical formula has only 1 enrichment value
			num_h, cf = NValueCalculator.parse_cf(partition[0])
			n_begin = 0
			low_pct = 0
			num_peaks = max(partition[1].intensities.str.len())
			emass_results_dict = dict()
			for row in partition[1].itertuples(index=False):
				enrichment = row.enrichment
				if enrichment == 0:
					enrichment = 0.05
				if enrichment not in emass_results_dict:
					high_pct = enrichment
					try:
						emass_results = emass(cf, n_begin, num_h, low_pct, high_pct, num_peaks)
					except Exception as e:
						print("Chemical Formula ", row.cf, " contains unsupported molecules")
					emass_results_dict[enrichment] = emass_results
				emass_results = emass_results_dict[enrichment]
				try:
					results.append(NValueCalculator.analyze_row(row, emass_results))
				except:
					print("EXCEPTION OCCURED WITH {}!".format(row))
			return pd.concat(
				[partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
				axis=1
			)
		except:
			return pd.concat(
				[partition[1].reset_index(drop=True), pd.DataFrame(data=results, columns=output_columns)],
				axis=1
			)

	@staticmethod
	def analyze_row(row, emass_results):
		"""
        Calculates n-value for each row

        Parameters:
            row (pd.Series): Contains chemical formula, empirical intensity data, and enrichment level,
            plot_dir (str): The directory to output fraction_new .csv files to.
            emass_results (pd.DataFrame): Results from emass containing unlabeled & labeled intensity and m/z data.

        Returns:
            Tuple: a tuple containing the given n-value and the stddev associated with it.
        """

		def s_n_filter(empirical_intensities, noise, NOISE_FILTER=10.0):
			empir_series = pd.Series(empirical_intensities)
			s_n_values = pd.Series([empirical_intensities[i] / noise for i in range(len(empirical_intensities))])
			filtered_empirical_intensities = empir_series.where(s_n_values >= NOISE_FILTER, np.nan)
			return filtered_empirical_intensities

		def dIt_filter(unfiltered_dataframe, dIt_data):
			return unfiltered_dataframe.where(dIt_data > 0.05, np.nan)

		def dIe_filter(unfiltered_dataframe, dIe_data):
			filtered_dataframe = unfiltered_dataframe.where(dIe_data < 0.05, np.nan)
			filtered_dataframe['n_D'] = unfiltered_dataframe['n_D']
			return filtered_dataframe

		def n_value_by_stddev(fraction_new, MINIMUM_PEAKS=3):
			stddev_dataframe = fraction_new.copy()

			# Find how many valid peaks exist
			stddev_dataframe = stddev_dataframe[
				['n_D'] + [col for col in stddev_dataframe.columns if 'I' in col]]
			stddev_dataframe['num_non_null'] = (stddev_dataframe.count(axis='columns') - 1)

			# rearrange columns to allow for proper std_dev calculation
			stddev_dataframe = stddev_dataframe[
				['num_non_null'] + [col for col in stddev_dataframe.columns if col != 'num_non_null']]

			# Calculate stddev for valid rows
			stddev_dataframe['stddev'] = stddev_dataframe.loc[:, 'I0'::].std(axis='columns')
			stddev_dataframe.loc[:, 'stddev'] = stddev_dataframe['stddev'].where(
				stddev_dataframe['num_non_null'] >= MINIMUM_PEAKS, np.nan)

			# Find n-value and stddev of filtered data
			if stddev_dataframe['stddev'].isna().all():
				filtered_nValue_min = np.nan
				filtered_stddev_min = np.nan

				truth_values = ~np.isnan(stddev_dataframe.iloc[:, 2:-1])
				truth_count = [np.sum(truth_values.iloc[:, x]) for x in range(truth_values.shape[1])]
				peaks_included = ' '.join(str(x) for x in truth_count)
			else:
				min_index = stddev_dataframe['stddev'].idxmin()
				filtered_nValue_min = stddev_dataframe.loc[min_index, 'n_D']
				filtered_stddev_min = stddev_dataframe['stddev'].min()

				data_mask = stddev_dataframe.iloc[min_index].notna().reset_index(drop=True)[2:-1]
				peaks_included = pd.Series(stddev_dataframe.columns)[2:-1][data_mask] \
					.to_string(header=False, index=False).replace('I', 'M').replace('\n', ' ')
				if peaks_included == 'Series([], )':
					peaks_included = ''
			return (stddev_dataframe, filtered_nValue_min, filtered_stddev_min, peaks_included)

		def n_value_dIt_filter(unfiltered_fraction_new, dIt_data, MINIMUM_PEAKS=3):
			###
			### Generate dIt filtered fraction_new DataFrame
			###
			filtered_fraction_new = dIt_filter(unfiltered_fraction_new, dIt_data)
			filtered_fraction_new.drop('stddev', axis=1, inplace=True)
			filtered_fraction_new.drop('num_non_null', axis=1, inplace=True)

			return n_value_by_stddev(filtered_fraction_new)

		def n_value_dIe_filter(unfiltered_fraction_new, dIe_data, MINIMUM_PEAKS=3):
			###
			### Generate dIt filtered fraction_new DataFrame
			###
			filtered_fraction_new = dIe_filter(unfiltered_fraction_new, dIe_data)
			filtered_fraction_new['n_D'] = unfiltered_fraction_new['n_D']
			filtered_fraction_new.drop('stddev', axis=1, inplace=True)
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
				'''Calculates the unit vector of a given vector
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
                '''

				def unit_vector(vector):
					'''Calculates the unit vector of a given vector
                    Parameters
                    ----------
                    vector : :obj:`list` of :obj:`float`
                    Returns
                    ----------
                    vector : :obj:`list` of :obj:`float`
                    '''
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

		# Returns 2 tuples of DataFrame with rows representing # of hydrogens that have been deuterated,
		# columns represent the peak #, data is m/z for 1st and intensity for 2nd
		(emass_unlabeled_mz, emass_unlabeled_intensities) = emass_unlabeled_data
		(emass_labeled_mz, emass_labeled_intensities) = emass_labeled_data

		unfiltered_fraction_new = emass_unlabeled_intensities.copy()

		# Generate fraction_new and dIt data
		for peak, ie in enumerate(row.intensities):
			unfiltered_fraction_new['I' + str(peak)] = (
					(emass_unlabeled_intensities['I' + str(peak)].iloc[0] - ie) /
					(emass_unlabeled_intensities['I' + str(peak)].iloc[0] - emass_labeled_intensities['I' + str(peak)])
			)

		# Calculate n-value with no filters (using stddev)
		unfiltered_fraction_new, n_value, stddev, peaks_included = n_value_by_stddev(unfiltered_fraction_new)

		return n_value, stddev

	@staticmethod
	def parse_cf(cf):
		d = dict(re.findall(r'([A-Z][a-z]*)(\d*)', cf))
		num_h = int(d['H'])
		blank = '{}'
		d.update({'H': blank, 'X': blank})  # H = # of Hydrogen, X = # labeled sites
		cf_string = ''.join('%s%s' % (k, v) for k, v in d.items())
		return num_h, cf_string


def main():
	print("This is a library file, please run the __main__ function in the same folder")


if __name__ == "__main__":
	main()
