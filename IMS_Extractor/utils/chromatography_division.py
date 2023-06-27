import pandas as pd
import numpy as np

# try:
from obs.peak import Peak
from obs.id import ID
from obs.envelope import Envelope
from obs.molecule import Molecule
import deuterater.settings as settings
# except:
# 	import IMS_Testing.obs.peak
# 	import IMS_Testing.obs.id
# 	import IMS_Testing.obs.envelope
# 	import IMS_Testing.obs.molecule
	# import IMS_Testing.deuterater.settings as settings
import multiprocessing as mp
import traceback
import os
from tqdm import tqdm
from functools import partial


class ChromatographyDivider:
	
	def __init__(self, settings_path, out_paths, input_paths, biomolecule_type):
		self.settings_path = settings_path
		settings.load(self.settings_path)
		
		self.input_paths = input_paths
		self.out_paths = out_paths
		self.how_divide = settings.use_chromatography_division
		self.biomolecule_type = biomolecule_type
		
		try:
			if settings.recognize_available_cores is True:
				self._n_processors = mp.cpu_count()
			else:
				self._n_processors = settings.n_processors
			if self._n_processors > 60:
				self._n_processors = 60
			self._mp_pool = mp.Pool(self._n_processors)
		
		except Exception as e:
			print(e)
			traceback.print_tb(e.__traceback__)
			raise
	
	@staticmethod
	def parse_2d_list(to_parse):
		split_list = to_parse.split("[")
		parsed = [split_list[i].replace("],", "").replace("]", "").split(", ")
				  for i in range(len(split_list))
				  if len(split_list[i]) > 1]
		
		parsed = [[parsed[i][j] for j in range(len(parsed[i])) if parsed[i][j] != ''] for i in range(len(parsed))]
		
		parsed = np.array(parsed).astype(float)
		
		return parsed
	
	@staticmethod
	def parse_1d_list(to_parse):
		parsed = to_parse.split(",")
		parsed[0] = parsed[0].replace("[", " ")
		parsed[-1] = parsed[-1].replace("]", "")
		
		return np.array(parsed)
	
	@staticmethod
	def divide_molecule(molecule):
		molecule.analyze_peaks(5)
		molecule.choose_peak()
		# molecule.extract_with_chosen_peak()
		return molecule
	
	@staticmethod
	def handle_molecule(group, settings_path):
		settings.load(settings_path)
		molecule_rt = Molecule(str(group[0]) + "_RT", settings_path)
		molecule_dt = Molecule(str(group[0]) + "_DT", settings_path)
		data = list()
		# Load in the data from the .tsv files
		for filename in np.unique(group[1]["infile"]):
			file_df = group[1].loc[group[1]["infile"] == filename]
			indexes = list(file_df["index"])
			indexes.append(-1)
			df = pd.read_csv(filename, sep="\t", skiprows=lambda x: x - 1 not in indexes)
			
			df["Extraction_Updated"] = ""
			df["Extraction_Error"] = ""
			df["mzml_name"] = df["mzml_path"].apply(os.path.basename)
			possible_rows = list(set(df.columns).intersection({"Adduct", "z", "mzml_name"}))
			possible_rows.sort()
			df['name_check'] = df.apply(lambda x: "_".join(x.loc[possible_rows].astype(str)), axis=1)
			df["filename"] = np.unique(file_df["filename"])[0]
			data.append(df)
		df = pd.concat(data)
		file_ids = dict()
		for series in df.iterrows():
			row = series[1]
			mzs = ChromatographyDivider.parse_2d_list(row.mzs_list)
			abs = ChromatographyDivider.parse_2d_list(row.intensities_list)
			rts = ChromatographyDivider.parse_1d_list(row.rt_list)
			dts = ChromatographyDivider.parse_1d_list(row.dt_list)
			baselines = ChromatographyDivider.parse_1d_list(row.baseline_list)
			
			# Get number of DT envelopes per RT:
			num_dt = len(np.unique(dts))
			num_rt = len(np.unique(rts))
			
			# Remove any data from memory so that I don't kill the memory.
			row.mzs_list = np.nan
			row.intensities_list = np.nan
			row.rt_list = np.nan
			row.dt_list = np.nan
			row.baseline_list = np.nan
			df.loc[df.name_check == row.name_check, "mzs_list"] = np.nan
			df.loc[df.name_check == row.name_check, "intensities_list"] = np.nan
			df.loc[df.name_check == row.name_check, "rt_list"] = np.nan
			df.loc[df.name_check == row.name_check, "dt_list"] = np.nan
			df.loc[df.name_check == row.name_check, "baseline_list"] = np.nan
			
			# Create envelopes that can be used to separate the data
			envelopes = list()
			should_plot = False
			num_normal_peaks = abs.shape[1] - settings.peak_lookback - settings.peak_lookahead
			for envelope_num in range(abs.shape[0]):
				try:
					peaks = [Peak(mzs[envelope_num][i], abs[envelope_num][i], i - 1) for i in range(abs.shape[1])]
				except:
					peaks = None
					print("hello")
				envelope = Envelope(peaks, rts[envelope_num], settings.peak_lookback, settings.peak_lookahead,
									dt=dts[envelope_num])
				envelope.baseline = baselines[envelope_num]
				
				def local_min(list_of_floats):
					for i in range(1, len(list_of_floats) - 1):
						if list_of_floats[i - 1] > list_of_floats[i] and list_of_floats[i + 1] > list_of_floats[i]:
							return False
					return True
				
				if len([peaks[i].ab for i in range(1, num_normal_peaks + 1)]) != num_normal_peaks or not local_min(
						[float(peak.ab) for peak in peaks[1:num_normal_peaks + 1]]):
					envelope.is_valid = False
				else:
					should_plot = True
				envelopes.append(envelope)
			
			id = ID(row.rt, row.mz, row.mass, row.z, row.n_isos)
			id.envelopes = envelopes
			envelope_2d = np.array(envelopes).reshape((num_rt, num_dt, len(envelopes[0])))
			
			# Separate a molecule by its Adduct, charge, and mzml_file
			column_list = row.index.isin(["Adduct", "z", "mzml_name"])
			column_list = row.index[column_list]
			column_list = list(column_list)
			column_list.sort()
			unique_identifier = "_".join(row.loc[column_list].astype(str))
			file_ids[unique_identifier] = id
			
			# print("Handle Dividing of RT and DT")
			all_abunds = np.vectorize((lambda x: x.ab))(envelope_2d)
			all_abunds = np.moveaxis(all_abunds, -1, 0)
			rts = rts.astype(float)
			dts = dts.astype(float)
			for i in range(num_dt):
				rt = rts[::num_dt]
				dt = dts[i]
				abunds = all_abunds[:, :, i]
				u_unique_identifier = f'{unique_identifier}_{dt}'
				molecule_rt.add_molecule_id(u_unique_identifier, rt, abunds, time_type="rt")

			for i in range(num_rt):
				rt = rts[::num_dt][i]
				dt = dts[:num_dt]
				abunds = all_abunds[:, i, :]
				u_unique_identifier = f'{unique_identifier}_{rt}'
				molecule_dt.add_molecule_id(u_unique_identifier, dt, abunds, time_type="dt")

		molecule_dt = ChromatographyDivider.divide_molecule(molecule_dt)
		molecule_rt = ChromatographyDivider.divide_molecule(molecule_rt)
		
		def sum_envelopes(envelopes):
			baselines = [e.baseline for e in envelopes]
			median_baseline = np.median(np.array(baselines).astype(float))
			rt = envelopes[0].rt
			peaks_ab = [[e._peaks[i].ab for e in envelopes] for i in range(len(envelopes[0]))]
			peaks_mz = [[e._peaks[i].mz for e in envelopes] for i in range(len(envelopes[0]))]
			summed_ab = np.array(peaks_ab).sum(axis=1)
			median_mz = np.median(np.array(peaks_mz), axis=1)
			peaks = [Peak(mz=median_mz[i], abundance=summed_ab[i], i=i-settings.peak_lookback) for i in range(len(summed_ab))]
			envelope = Envelope(peaks, rt=rt, n_lookback=settings.peak_lookback, n_lookahead=settings.peak_lookahead)
			envelope.baseline = median_baseline
			return envelope

		def update_output_file(df, id, name, dt_min, dt_max):
			df["row_num"] = np.arange(0, df.shape[0])
			df.set_index("row_num", drop=False, inplace=True)
			df['signal_2_noise'] = ''
			df['signal_2_noise'].astype(str)
			df["summed_signal"] = df["summed_signal"].astype(str)
			to_string = ['m-1_mz', 'm-1_abundance', 'm_end+1_mz', 'm_end+1_mz', 'm_end+1_abundance', 'mads', 'mzs',
						 'abundances']
			df = df.astype({a: str for a in to_string})
			name_column = list(set(df.columns).intersection({"Lipid Unique Identifier", "Sequence"}))
			# row = df[df[name_column[0]] == name]
			row = df.loc[df["name_check"] == name]
			for row_index in row["row_num"]:
				if id.condensed_envelope:
					mzs, abundances = id.condensed_envelope.to_obs()
					summed_signal = str(id.total_signal)
					signal_2_noise = [float(a) / float(id.condensed_envelope.baseline) for a in abundances]
					lb_mzs, lb_abundances = id.condensed_envelope.lb_obs()
					la_mzs, la_abundances = id.condensed_envelope.la_obs()
					df.at[row_index, 'mzs'] = str(mzs)
					df.at[row_index, 'abundances'] = str(abundances)
					df.at[row_index, 'summed_signal'] = summed_signal
					df.at[row_index, 'rt_min'] = id.rt_min
					df.at[row_index, 'rt_max'] = id.rt_max
					df.at[row_index, "Extraction_Updated"] = "OUCH"  # self.why_chosen
					df.at[row_index, 'baseline_signal'] = id.condensed_envelope.baseline
					df.at[row_index, 'signal_2_noise'] = str(signal_2_noise)
					df.at[row_index, 'm-1_mz'] = str(lb_mzs)
					df.at[row_index, 'm-1_abundance'] = str(lb_abundances)
					df.at[row_index, 'm_end+1_mz'] = str(la_mzs)
					df.at[row_index, 'm_end+1_abundance'] = str(la_abundances)
					df.at[row_index, 'mads'] = str(id.mads)
					df.at[row_index, 'num_scans_combined'] = len(id.envelopes)
					df.at[row_index, 'dt_min'] = dt_min
					df.at[row_index, 'dt_max'] = dt_max
					df.at[row_index, 'updated'] = True
				else:
					try:
						df.at[row_index, 'mzs'] = np.nan
					except:
						print("There is an errror in the else")
						continue
					df.at[row_index, 'abundances'] = np.nan
					df.at[row_index, 'summed_signal'] = "[]"
					df.at[row_index, 'rt_min'] = np.nan
					df.at[row_index, 'rt_max'] = np.nan
					df.at[row_index, 'dt_min'] = np.nan
					df.at[row_index, 'dt_max'] = np.nan
					df.at[row_index, 'baseline_signal'] = np.nan
					df.at[row_index, 'signal_2_noise'] = np.nan
					df.at[row_index, 'm-1_mz'] = np.nan
					df.at[row_index, 'm-1_abundance'] = np.nan
					df.at[row_index, 'm_end+1_mz'] = np.nan
					df.at[row_index, 'm_end+1_abundance'] = np.nan
					df.at[row_index, 'mads'] = np.nan
					df.at[row_index, 'num_scans_combined'] = 0
					df.at[row_index, "Extraction_Updated"] = "no peak used"
					df.at[row_index, "Extraction_Error"] = "An Error Occurred in Chromatography Division"
					df.at[row_index, 'updated'] = True
					# self.error[name]
			return df

		try:
			dt_indexes = molecule_dt.unique_chrom_peaks[molecule_dt.chosen_peak_index]
			dt_min = dts[dt_indexes[0]]
			dt_max = dts[dt_indexes[1]]
			rt_indexes = molecule_rt.unique_chrom_peaks[molecule_rt.chosen_peak_index]
		except:
			dt_indexes = (0, 0)
			dt_min = None
			dt_max = None
			rt_indexes = (0, 0)

		updated_dfs = list()
		
		for name, id in file_ids.items():
			if dt_min is None:
				df = update_output_file(df, id, name, dt_min, dt_max)
				continue
			envelopes = id.envelopes
			summed_envelopes = list()
			for i in range(*rt_indexes):
				used_envelopes = envelopes[num_dt*i+dt_indexes[0]:num_dt*i+dt_indexes[1]]
				summed_envelopes.append(sum_envelopes(used_envelopes))
			id.envelopes = summed_envelopes
			id.aggregate_envelopes()
			df = update_output_file(df, id, name, dt_min, dt_max)
		group_df = df
		# group_df = pd.concat(updated_dfs)
		columns_to_drop = ["name_check", "mzml_name", "row_num",
						   "mzs_list", "intensities_list", "rt_list", "baseline_list"]
		column_renames = {"Extraction_Error": "Chromatography_Division_Error",
						  "Extraction_Updated": "Chromatography_Division_Scores", }
		
		columns_to_use = [x for x in list(group_df.columns) if x not in columns_to_drop]
		group_df = group_df[columns_to_use].rename(columns=column_renames)
		
		return group_df
	
	def divide(self):
		if self.biomolecule_type == "Lipid":
			col_names = ["Lipid Unique Identifier", 'Precursor Drift Time (ms)', "rt_list", "mzml_path", "Adduct", "z"]
			molecule_group_name = "Lipid Unique Identifier"
		elif self.biomolecule_type == "Peptide":
			col_names = ["Sequence", "rt_list", "mzml_path", "z"]
			molecule_group_name = "Sequence"
		else:
			col_names = []
			molecule_group_name = ""
		if self.how_divide == "Interfile":
			dataframes = list()
			for i in range(len(self.input_paths)):
				infile, outfile = self.input_paths[i], self.out_paths[i]
				df = pd.read_csv(infile, sep="\t", usecols=col_names)
				df = df.dropna(subset=["rt_list"])
				del df["rt_list"]
				df["filename"] = outfile
				df["infile"] = infile
				dataframes.append(df)
			df = pd.concat(dataframes)
			df = df.reset_index()
			df["mzml_name"] = df["mzml_path"].apply(os.path.basename)
			df[molecule_group_name] = df[molecule_group_name] + df["Precursor Drift Time (ms)"].astype(str)
			molecule_groups = df.groupby(by=molecule_group_name)
			molecules = list()
			
			num_groups = molecule_groups.ngroups
			if settings.debug_level == 0:
				func = partial(self.handle_molecule, settings_path=self.settings_path)
				molecules = self._mp_pool.map(func,
											  tqdm(molecule_groups, desc="dividing chromatography: ",
												   total=num_groups))
			
			elif settings.debug_level >= 1:
				func = partial(self.handle_molecule, settings_path=self.settings_path)
				cool = True
				for group in tqdm(molecule_groups, desc="dividing chromatography: ", total=num_groups):
					molecules.append(func(group))
			
			df = pd.concat(molecules)
			
			file_groups = df.groupby(by="filename")
			for file_group_name, file_group in file_groups:
				file_group.to_csv(file_group_name, sep="\t", index=False)
		
		else:
			for i in tqdm(range(len(self.out_paths)), total=len(self.out_paths), desc="dividing files: "):
				infile, outfile = self.input_paths[i], self.out_paths[i]
				df = pd.read_csv(infile, sep='\t', usecols=col_names)
				df['filename'] = outfile
				df["infile"] = infile
				df = df.reset_index()
				df = df.dropna(subset=["rt_list"])
				del df["rt_list"]
				# possible_rows = list(set(df.columns).intersection({"Adduct", "z", "mzml_name"}))
				# possible_rows.sort()
				# df['name_check'] = df.apply(lambda x: "_".join(x.loc[possible_rows].astype(str)), axis=1)
				molecule_groups = df.groupby(by=molecule_group_name)
				# df["Extraction_Updated"] = ""
				# df["Extraction_Error"] = ""
				molecules = list()
				
				if settings.debug_level == 0:
					func = partial(self.handle_molecule, settings_path=self.settings_path)
					molecules = self._mp_pool.map(func,
												  tqdm(molecule_groups,
													   desc="dividing chromatography: ",
													   total=len(molecule_groups),
													   leave=False))
				
				elif settings.debug_level >= 1:
					func = partial(self.handle_molecule, settings_path=self.settings_path)
					for group in tqdm(molecule_groups, desc="dividing chromatography: ", total=len(molecule_groups),
									  leave=False):
						molecules.append(func(group))
				
				df = pd.concat(molecules)
				
				df.to_csv(outfile, sep='\t', index=False)


def main():
	input_files = [
		"/media/JCPriceLab/BUBBLES/Chad Q/DR++IMS/IMS_Extractor/test1_only_extracted_sorted32_5dtw.tsv"]
	out_files = [a[:-4] + "_divided.tsv" for a in input_files]
	settings_file = "/media/JCPriceLab/BUBBLES/Chad Q/DR++IMS/IMS_Extractor/settings.yaml"
	# df = pd.read_csv("/home/JCPriceLab/Documents/PriceLabSoftware/DeuteRater Extensions/IMS-Testing/D0_MS_only_extracted.tsv", sep='\t')
	divider = ChromatographyDivider(settings_path=settings_file,
									input_paths=input_files,
									out_paths=out_files,
									biomolecule_type="Lipid")
	divider.divide()


# print("This is a utility file. ",
# 	  "If you would like to use the Chromatography Division by itself, ",
# 	  "run the file chromatography_division_sep.py inside the main DeuteRater folder.")


if __name__ == "__main__":
	main()
