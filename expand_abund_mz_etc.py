import pandas as pd
import numpy as np
from tqdm import tqdm

def tk_get_files(extension='*', prompt="Select file"):
	from tkinter import filedialog
	from tkinter import Tk
	import os
	root = Tk()
	root.withdraw()
	if (extension == '*'):
		root.filename = filedialog.askopenfilenames(
			initialdir=os.getcwd(), title=prompt)
	else:
		extension_list = list()
		extension.split(",")
		for extension in extension.split(","):
			if extension == ' ':
				continue
			elif extension[0] == ' ':
				extension = extension[1:]
			elif extension[0] != '.':
				extension = "." + extension
			extension_list.append((extension + " Files", extension), )
		extension_list.append(("All Files", "*"), )
		
		root.filename = filedialog.askopenfilenames(
			initialdir=os.getcwd(), title=prompt,
			filetypes=extension_list)
	root.update()
	
	filename = root.filename
	if len(filename) == 0:
		print("No File Grabbed")
		exit(1)
	return filename  # A string representing the file path

fix_signal_2_noise = False

# filename = "/mnt/wwn-0x5000c500c4ea04a1-part2/CQ1/CQ1105/Corrected mz/A2_pos/Change Column Names 18 Oct/A2D0M1_pos.tsv"
filenames = tk_get_files()
for filename in tqdm(filenames):
	
	df = pd.read_csv(filename, sep="\t")
	
	to_sep = ["mzs", "abundances", "signal_2_noise", "mads"]
	
	fn_seps = ["theory_unlabeled_abunds", "theory_labeled_abunds", "normalized_empirical_abundances", "frac_new_abunds"]
	
	if fix_signal_2_noise:
		try:
			df["baseline_signal"] = df["baseline_signal"] * 3
			a = df["abundances"]
			# Turn abundances to a list of values:
			a = a.dropna().apply(lambda x: np.array(x[1:-2].split(", ")).astype(float))
			a = a / df["baseline_signal"].dropna()
			df["signal_2_noise"] = a
		except Exception as e:
			print("singal 2 noise error")
		
	if "frac_new_abunds_std_dev" in df.columns:
		to_sep = to_sep + fn_seps
	
	for col_name in to_sep:
		# if col_name == "signal_2_noise":
		# 	separated_cols = pd.DataFrame(df[col_name]
		# 								  .tolist(),
		# 								  columns=[f"{col_name}_m0", f"{col_name}_m1", f"{col_name}_m2"])
		# 	df = df.join(separated_cols)
		# else:
			try:
				separated_cols = pd.DataFrame(df[col_name]
											  .apply(lambda x: str(x)[1:-1].split(",") if str(x) != "nan" else ["nan", "nan", "nan"])
											  .tolist(),
											  columns=[f"{col_name}_m0", f"{col_name}_m1", f"{col_name}_m2"])
				df = df.join(separated_cols)
			except Exception as e:
				print("expansion error")
	
	if "Extraction Updated" in df.columns:
		try:
			scores = df["Extraction Updated"].split(", ")
			df["Adduct Score"] = scores[0]
			df["Neutromer Peak Tracing Score"] = scores[1]
			df["Distance from Reported RT Score"] = scores[2]
			df["Intensity Score"] = scores[3]
		except:
			pass
	
	to_sep_hi_lo = ["dietary_lit_n_range", "de_novo_lit_n_range"]
	for col_name in to_sep_hi_lo:
		try:
			separated_cols = pd.DataFrame(df[col_name]
										  .apply(
				lambda x: str(x)[1:-1].split(", ") if str(x) != "nan" else ["nan", "nan"])
										  .tolist(),
										  columns=[f"{col_name}_low", f"{col_name}_hi"])
		except:
			pass
	
	df.rename({"baseline_signal": "noise"}, inplace=True)
	
	df.to_csv(filename[:-4] + "_splitup.tsv", sep="\t", index=False)
