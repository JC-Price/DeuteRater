import os
from functools import partial

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# biomolecule: str


def tk_get_single_file(extension='*', prompt="Select file"):
	from tkinter import filedialog
	from tkinter import Tk
	import os
	if extension[0] != '*':
		extension = '*' + extension
	root = Tk()
	root.withdraw()
	if (extension == '*'):
		root.filename = filedialog.askopenfilenames(
			initialdir=os.getcwd(), title=prompt)
	else:
		root.filename = filedialog.askopenfilenames(
			initialdir=os.getcwd(), title=prompt,
			filetypes=((extension[2:] + " files", extension)))
	root.update()
	
	filename = root.filename
	return filename  # A string representing the file path


def tk_get_all_in_folder(extension='*', prompt="Select folder"):
	from tkinter import filedialog
	from tkinter import Tk
	import fnmatch
	import os
	if extension[0] != '*':
		extension = '*' + extension
	root = Tk()
	root.withdraw()
	root.filename = filedialog.askdirectory(
		initialdir=os.getcwd(), title=prompt)
	root.update()
	folder = root.filename
	
	all_folder_files = os.listdir(folder)
	
	path_sep = os.path.sep
	files_in_folder = list()
	for file in all_folder_files:
		if fnmatch.fnmatch(file, extension):
			files_in_folder.append(folder + path_sep + file)
	
	return files_in_folder  # a list of file paths that have the required extension


# def get_file_path_from_file(file):
# 	import os
# 	return os.path.sep.join(file.split(os.path.sep)[:-1]) + os.path.sep


# def get_file_name_from_path(file):
# 	import os
# 	return file.split(os.path.sep)[-1]


def parse_2d_list(to_parse):
	split_list = to_parse.split("[")
	parsed = [split_list[i].replace("]", "").split(", ")
			  for i in range(len(split_list))
			  if len(split_list[i]) > 1]
	
	parsed = [[parsed[i][j] for j in range(len(parsed[i])) if parsed[i][j] != ''] for i in range(len(parsed))]
	
	parsed = np.array(parsed).astype(float)
	
	return parsed


def parse_1d_list(to_parse):
	parsed = to_parse.split(",")
	parsed[0] = parsed[0].replace("[", "")
	parsed[-1] = parsed[-1].replace("]", "")
	
	return np.array(parsed)


def plot_subplot(axs, data, intensity_filter=20_000, rt_list=None, nisos=3, Adduct=""):
	def format_for_empty(axs, rt_list, Adduct, nisos):
		rt_list = parse_1d_list(rt_list)
		m0 = np.zeros(rt_list.shape)
		m0[0] = 100_000_000
		combined = m0
		mz_value = "no data available"
		axs.plot(rt_list, m0, label="m0")
		
		for i in range(1, nisos):
			mi = np.zeros(rt_list.shape)
			axs.plot(rt_list, mi, label=f"m{i}")
		
		axs.plot(rt_list, combined, label="combined")
		axs.set_title(str(Adduct) + "\n m/z = " + mz_value)
		axs.grid(axis="y")
		axs.set_yticks(np.linspace(0, np.max(combined), 9))
		axs.axvline(x=0, c='k', label="Extracted RT Range")
		axs.axvline(x=0, c='k')
		axs.set_facecolor("darkslategrey")
		return axs
	
	if data.empty:
		return format_for_empty(axs, rt_list, Adduct, nisos)
	
	row = data.iloc[0]
	# mz_list = parse_2d_list(row.mzs_list)
	ab_list = parse_2d_list(row.intensities_list)
	rt_list = parse_1d_list(row.rt_list)
	
	m_lists = list()
	for m_val in range(1, ab_list.shape[1] - 1):
		m_lists.append(ab_list[:, m_val])
	
	# m_1 = ab_list[:, 0]
	max_intensities = [max(m) for m in m_lists]
	
	combined = sum(m_lists)
	
	if combined.max == 0.0:
		return format_for_empty(axs, rt_list, Adduct)
	
	for m_val in range(nisos):
		axs.plot(rt_list, m_lists[m_val], label=f"m{m_val}")
	axs.plot(rt_list, combined, label="combined")
	# axs.plot(rt_list, m_1, label="m-1")
	
	axs.set_title(str(Adduct) + "\n m/z = " + str(row["Precursor m/z"]))
	axs.grid(axis="y")
	
	# Set the yticks so every plot has 5
	axs.set_yticks(np.linspace(0, np.max(combined), 9))
	
	# Put lines around rt_start and rt_end if they exist
	def rt_index(rt_list, value):
		return (np.abs(rt_list.astype(float) - float(value))).argmin()
	
	if not np.isnan(row.rt_min):
		axs.axvline(x=rt_index(rt_list, row.rt_min), c='k', label="Extracted RT Range")
		axs.axvline(x=rt_index(rt_list, row.rt_max), c='k')
	else:
		axs.axvline(x=0, c='k', label="Extracted RT Range")
		axs.axvline(x=0, c='k')
		if max(max_intensities) < intensity_filter:
			axs.set_facecolor("powderblue")
		elif "Extraction_Error" in data.columns:
			if data.iloc[0]["Extraction_Error"] == "Chromatographic Peak does not exist for this rep.":
				axs.set_facecolor("lightgray")
			elif data.iloc[0]["Extraction_Error"] == "CFS, Angle, Threshold Filtered resulted in less than 10 scans.":
				axs.set_facecolor("yellow")
		else:
			axs.set_facecolor("pink")
	
	return axs


def handle_single_file(lipids_selected, biomolecule_type, files):
	def parse_2d_list(to_parse):
		split_list = to_parse.split("[")
		parsed = [split_list[i].replace("]", "").split(", ")
				  for i in range(len(split_list))
				  if len(split_list[i]) > 1]
		
		parsed = [[parsed[i][j] for j in range(len(parsed[i])) if parsed[i][j] != ''] for i in range(len(parsed))]
		
		parsed = np.array(parsed).astype(float)
		
		return parsed
	
	def parse_1d_list(to_parse):
		parsed = to_parse.split(",")
		parsed[0] = parsed[0].replace("[", "")
		parsed[-1] = parsed[-1].replace("]", "")
		
		return np.array(parsed)
	
	def plot_subplot(axs, data, intensity_filter=20_000, rt_list=None, nisos=3, Adduct=""):
		def format_for_empty(axs, rt_list, Adduct, nisos):
			rt_list = parse_1d_list(rt_list)
			m0 = np.zeros(rt_list.shape)
			m0[0] = 100_000_000
			combined = m0
			mz_value = "no data available"
			axs.plot(rt_list, m0, label="m0")
			
			for i in range(1, nisos):
				mi = np.zeros(rt_list.shape)
				axs.plot(rt_list, mi, label=f"m{i}")
			
			axs.plot(rt_list, combined, label="combined")
			axs.set_title(str(Adduct) + "\n m/z = " + mz_value)
			axs.grid(axis="y")
			axs.set_yticks(np.linspace(0, np.max(combined), 9))
			axs.axvline(x=0, c='k', label="Extracted RT Range")
			axs.axvline(x=0, c='k')
			axs.set_facecolor("darkslategrey")
			return axs
		
		if data.empty:
			return format_for_empty(axs, rt_list, Adduct, nisos)
		
		row = data.iloc[0]
		# mz_list = parse_2d_list(row.mzs_list)
		ab_list = parse_2d_list(row.intensities_list)
		rt_list = parse_1d_list(row.rt_list)
		
		m_lists = list()
		for m_val in range(1, ab_list.shape[1] - 1):
			m_lists.append(ab_list[:, m_val])
		
		# m_1 = ab_list[:, 0]
		max_intensities = [max(m) for m in m_lists]
		
		combined = sum(m_lists)
		
		if combined.max == 0.0:
			return format_for_empty(axs, rt_list, Adduct)
		
		for m_val in range(nisos):
			axs.plot(rt_list, m_lists[m_val], label=f"m{m_val}")
		axs.plot(rt_list, combined, label="combined")
		# axs.plot(rt_list, m_1, label="m-1")
		
		axs.set_title(str(Adduct) + "\n m/z = " + str(row["Precursor m/z"]))
		axs.grid(axis="y")
		
		# Set the yticks so every plot has 5
		axs.set_yticks(np.linspace(0, np.max(combined), 9))
		
		# Put lines around rt_start and rt_end if they exist
		def rt_index(rt_list, value):
			return (np.abs(rt_list.astype(float) - float(value))).argmin()
		
		if not np.isnan(row.rt_min):
			axs.axvline(x=rt_index(rt_list, row.rt_min), c='k', label="Extracted RT Range")
			axs.axvline(x=rt_index(rt_list, row.rt_max), c='k')
		else:
			axs.axvline(x=0, c='k', label="Extracted RT Range")
			axs.axvline(x=0, c='k')
			if max(max_intensities) < intensity_filter:
				axs.set_facecolor("powderblue")
			elif "Extraction_Error" in data.columns:
				if data.iloc[0]["Extraction_Error"] == "Chromatographic Peak does not exist for this rep.":
					axs.set_facecolor("lightgray")
				elif data.iloc[0][
					"Extraction_Error"] == "CFS, Angle, Threshold Filtered resulted in less than 10 scans.":
					axs.set_facecolor("yellow")
			else:
				axs.set_facecolor("pink")
		
		return axs
	
	# global biomolecule_type
	
	chrom_file = files[0]
	divided_file = files[1]
	if chrom_file[-4:] != ".tsv" and divided_file[-4:] != ".tsv":
		return
	
	# Before reading in the chromatography data, read in all available lipids to find the
	# index that the needed lipids are located in the file being read in to save time and space.
	if biomolecule_type == "Lipid":
		col_type = "Lipid Unique Identifier"
		index_col = "Adduct_cf"
	else:
		col_type = "Sequence"
		index_col = "cf"
	
	available_lipids = pd.read_csv(divided_file, usecols=[col_type], sep='\t')
	
	for lipid_selected in lipids_selected:
		folder_name = lipid_selected
		folder_name = folder_name.replace(":", "_")
		folder_name = folder_name.replace("/", "_")
		
		lipid_name = lipid_selected
		lipid_name = lipid_name.replace(":", "_")
		lipid_name = lipid_name.replace("/", "_")
		
		# Create a folder to store the images:
		folder_path = os.path.join(os.path.dirname(files[0]), folder_name)
		try:
			if not os.path.isdir(folder_path):
				os.mkdir(folder_path)
		except OSError as e:
			pass
		
		rt_indexes = available_lipids.loc[available_lipids[col_type] == lipid_selected].index
		if len(rt_indexes) == 0:
			print("That lipid or charge does not exist.")
			continue
		
		rt_indexes = [index + 1 for index in rt_indexes]
		rt_indexes.insert(0, 0)  # This ensures the heading is read in correctly.
		
		if divided_file == chrom_file:
			chrom_lipids = available_lipids.copy()
		else:
			chrom_lipids = pd.read_csv(chrom_file, usecols=[col_type], sep='\t')
		lipid_indexes = chrom_lipids.loc[chrom_lipids[col_type] == lipid_selected].index
		lipid_indexes = [index + 1 for index in lipid_indexes]
		lipid_indexes.insert(0, 0)
		
		# Read in just the rows with data we need.
		lipid = pd.read_csv(chrom_file, sep='\t', skiprows=(lambda x: x not in lipid_indexes))
		if chrom_file == divided_file and rt_indexes == lipid_indexes:
			rt_data = lipid.copy()
		else:
			rt_data = pd.read_csv(divided_file, sep='\t', skiprows=(lambda x: x not in rt_indexes))
		lipid = lipid.set_index(index_col)
		rt_data = rt_data.set_index(index_col)
		
		# grab the RT from the divided file
		lipid['rt_min'] = rt_data['rt_min']
		lipid['rt_max'] = rt_data['rt_max']
		
		lipid = lipid.loc[~lipid["rt_list"].isna()]  # Remove any rows where the lipid was not actually extracted.
		if lipid.empty:
			plt.title(f"{lipid_selected} was never extracted")
			plt.savefig(os.path.join(folder_path, os.path.basename(divided_file)[:-4] + "_" + lipid_name + ".png"))
			plt.close()
			continue
		
		if biomolecule_type == "Lipid":
			fig, axs = plt.subplots(6, sharex=True, figsize=(10, 22))
		else:
			fig, axs = plt.subplots(5, sharex=True, figsize=(10, 22))
		
		adduct_order = ["M+H", "M+H-[H20]", "M+Na",
						"M+Na-[H2O]", "M+NH4", "M+NH4-[H2O]"]
		
		charge_order = [1, 2, 3, 4, 5]
		
		rt_list = lipid["rt_list"].iloc[0]
		
		nisos = lipid["n_isos"].iloc[0]
		
		if biomolecule_type == "Lipid":
			for i in range(len(adduct_order)):
				axs[i] = plot_subplot(axs[i], lipid.loc[lipid["Adduct"] == adduct_order[i]], rt_list=rt_list,
									  nisos=nisos, Adduct=adduct_order[i])
		else:
			for i in range(len(charge_order)):
				axs[i] = plot_subplot(axs[i], lipid.loc[lipid["z"] == charge_order[i]], rt_list=rt_list, nisos=nisos,
									  Adduct=str(charge_order[i]))
		
		rt_list = parse_1d_list(lipid.iloc[0].rt_list).astype(float)
		
		fig_title = f"{lipid_selected} z={1}"
		if "Extraction_Updated" in lipid.columns and len(np.unique(lipid["Extraction_Updated"])) > 1:
			fig_title = f"{fig_title}\n{[a for a in np.unique(lipid['Extraction_Updated']) if a != 'no peak used'][0]}"
		fig.suptitle(fig_title, fontsize=20)
		
		plt.xticks(range(0, len(rt_list), int(len(rt_list) / 180)), rt_list[0::int(len(rt_list) / 180)], rotation=90)
		handles, labels = axs[-1].get_legend_handles_labels()
		
		from matplotlib.patches import Patch
		pink_patch = Patch(color="pink", label="other issues")
		blue_patch = Patch(color="powderblue", label="too low of intensity")
		yellow_patch = Patch(color="yellow", label="CFS/Angle Issue")
		gray_patch = Patch(color="lightgray", label="No Peak Available")
		no_data_patch = Patch(color="darkslategrey", label="No Data Extracted")
		
		handles.append(blue_patch)
		labels.append(blue_patch.get_label())
		handles.append(yellow_patch)
		labels.append(yellow_patch.get_label())
		handles.append(gray_patch)
		labels.append(gray_patch.get_label())
		handles.append(pink_patch)
		labels.append(pink_patch.get_label())
		handles.append(no_data_patch)
		labels.append(no_data_patch.get_label())
		
		fig.legend(handles, labels, ncol=4, fontsize=15, loc="upper center", bbox_to_anchor=(0.5, 0.95))
		
		# plt.ylabel("Intensity")
		# plt.xlabel("RT")
		
		# plt.show()
		
		plt.savefig(os.path.join(folder_path, os.path.basename(divided_file)[:-4] + "_" + lipid_name + ".png"))
		plt.close()


def handle_single_lipid(lipid_selected, chrom_files, divided_files):
	folder_name = lipid_selected
	folder_name = folder_name.replace(":", "_")
	folder_name = folder_name.replace("/", "_")
	
	lipid_name = lipid_selected
	lipid_name = lipid_name.replace(":", "_")
	lipid_name = lipid_name.replace("/", "_")
	
	for i in range(len(chrom_files)):
		chrom_file = chrom_files[i]
		divided_file = divided_files[i]
		if chrom_file[-4:] != ".tsv" and divided_file[-4:] != ".tsv":
			continue
		# Create a folder to store the images:
		folder_path = os.path.dirname(divided_file) + folder_name + os.path.sep
		try:
			os.mkdir(folder_path)
		except OSError as e:
			pass
		
		# Before reading in the chromatography data, read in all available lipids to find the
		# index that the needed lipids are located in the file being read in to save time and space.
		global biomolecule_type
		if biomolecule_type == "Lipid":
			col_type = "Lipid Unique Identifier"
			index_col = "Adduct_cf"
		else:
			col_type = "Sequence"
			index_col = "cf"
		
		available_lipids = pd.read_csv(divided_file, usecols=[col_type], sep='\t')
		rt_indexes = available_lipids.loc[available_lipids[col_type] == lipid_selected].index
		if len(rt_indexes) == 0:
			print("That lipid or charge does not exist.")
			continue
		rt_indexes = [index + 1 for index in rt_indexes]
		rt_indexes.insert(0, 0)  # This ensures the heading is read in correctly.
		
		chrom_lipids = pd.read_csv(chrom_file, usecols=[col_type], sep='\t')
		lipid_indexes = chrom_lipids.loc[chrom_lipids[col_type] == lipid_selected].index
		lipid_indexes = [index + 1 for index in lipid_indexes]
		lipid_indexes.insert(0, 0)
		
		# Read in just the rows with data we need.
		lipid = pd.read_csv(chrom_file, sep='\t', skiprows=(lambda x: x not in lipid_indexes))
		rt_data = pd.read_csv(divided_file, sep='\t', skiprows=(lambda x: x not in rt_indexes))
		lipid = lipid.set_index(index_col)
		rt_data = rt_data.set_index(index_col)
		
		# grab the RT from the divided file
		lipid['rt_min'] = rt_data['rt_min']
		lipid['rt_max'] = rt_data['rt_max']
		
		lipid = lipid.loc[~lipid["rt_list"].isna()]  # Remove any rows where the lipid was not actually extracted.
		if lipid.empty:
			plt.title(f"{lipid_selected} was never extracted")
			plt.savefig(folder_path + os.path.basename(divided_file)[:-4] + "_" + lipid_name + ".png")
			plt.close()
			continue
		
		if biomolecule_type == "Lipid":
			fig, axs = plt.subplots(6, sharex=True, figsize=(10, 22))
		else:
			fig, axs = plt.subplots(5, sharex=True, figsize=(10, 22))
		
		adduct_order = ["M+H", "M+H-[H20]", "M+Na",
						"M+Na-[H2O]", "M+NH4", "M+NH4-[H2O]"]
		
		charge_order = [1, 2, 3, 4, 5]
		
		rt_list = lipid["rt_list"].iloc[0]
		
		nisos = lipid["n_isos"].iloc[0]
		
		if biomolecule_type == "Lipid":
			for i in range(len(adduct_order)):
				axs[i] = plot_subplot(axs[i], lipid.loc[lipid["Adduct"] == adduct_order[i]], rt_list=rt_list,
									  nisos=nisos, Adduct=adduct_order[i])
		else:
			for i in range(len(charge_order)):
				axs[i] = plot_subplot(axs[i], lipid.loc[lipid["z"] == charge_order[i]], rt_list=rt_list, nisos=nisos,
									  Adduct=str(charge_order[i]))
		
		rt_list = parse_1d_list(lipid.iloc[0].rt_list).astype(float)
		
		fig_title = f"{lipid_selected} z={1}"
		if "Extraction_Updated" in lipid.columns and len(np.unique(lipid["Extraction_Updated"])) > 1:
			fig_title = f"{fig_title}\n{[a for a in np.unique(lipid['Extraction_Updated']) if a != 'no peak used'][0]}"
		fig.suptitle(fig_title, fontsize=20)
		# axs.set_xticks(np.arange(len(rt_list), 10))
		# 	axs.set_xticklabels(rt_list[0::10])
		plt.xticks(range(0, len(rt_list), 10), rt_list[0::10], rotation=90)
		handles, labels = axs[-1].get_legend_handles_labels()
		
		from matplotlib.patches import Patch
		pink_patch = Patch(color="pink", label="other issues")
		blue_patch = Patch(color="powderblue", label="too low of intensity")
		yellow_patch = Patch(color="yellow", label="CFS/Angle Issue")
		gray_patch = Patch(color="lightgray", label="No Peak Available")
		no_data_patch = Patch(color="darkslategrey", label="No Data Extracted")
		
		handles.append(blue_patch)
		labels.append(blue_patch.get_label())
		handles.append(yellow_patch)
		labels.append(yellow_patch.get_label())
		handles.append(gray_patch)
		labels.append(gray_patch.get_label())
		handles.append(pink_patch)
		labels.append(pink_patch.get_label())
		handles.append(no_data_patch)
		labels.append(no_data_patch.get_label())
		
		fig.legend(handles, labels, ncol=4, fontsize=15, loc="upper center", bbox_to_anchor=(0.5, 0.95))
		
		# plt.ylabel("Intensity")
		# plt.xlabel("RT")
		
		# plt.show()
		
		plt.savefig(folder_path + os.path.basename(divided_file)[:-4] + "_" + lipid_name + ".png")
		plt.close()


def main():
	# global biomolecule_type
	biomolecule_type = input("What biomolecule did you use? (Lipid/Peptide): ")
	
	if biomolecule_type[0] == "l" or biomolecule_type == "L":
		biomolecule_type = "Lipid"
	else:
		biomolecule_type = "Peptide"
	
	lipids_selected = input(
		f"Please input all {biomolecule_type}s Unique Identifiers you want, separated by a comma and space ', ': ")
	lipids = lipids_selected.split(", ")
	charge_selected = 1
	
	division_done = input("Was chromatography division performed (y/n): ")
	
	if division_done == "y" or division_done == "Y":
		chromatography_files = list(tk_get_single_file(prompt="Select no_division Files: "))
		divided_files = list(tk_get_single_file(prompt="Select Corresponding Divided Files: "))
		chromatography_files.sort()
		divided_files.sort()
	else:
		chromatography_files = tk_get_single_file(prompt="Select Extracted Files: ")
		divided_files = chromatography_files
	
	# for lipid in lipids:
	# 	handle_single_lipid(lipid, chromatography_files, divided_files)
	#
	from tqdm import tqdm
	import multiprocessing as mp
	
	try:
		cpus = mp.cpu_count()
	except:
		cpus = 2  # Default value
	files = [(chromatography_files[i], divided_files[i]) for i in range(len(chromatography_files))]
	for file in tqdm(files):
		handle_single_file(lipids, biomolecule_type=biomolecule_type, files=file)
	# with mp.Pool(processes=cpus) as pool:
	# 	func = partial(handle_single_file, lipids_selected=lipids)
	# 	func = partial(func, biomolecule_type=biomolecule_type)
	# 	pool.map(func, tqdm(files))


if __name__ == "__main__":
	main()
