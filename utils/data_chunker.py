from utils.NValueCalculator import NValueCalculator
import time
import pandas as pd
import string

from tkinter import filedialog
from tkinter import Tk
import os

output_columns = ['n_value', 'stddev']

def data_chunker(df_full):
	# Renames columns to proper names depending on if a lipid or protein is being passed in
	df = df_full.rename(columns={  # TODO: Verify Column Names
		'cf': 'chemical_formula',
		'abundances': 'intensities',
		'labeling': 'enrichment'})

	df['intensities'] = df['intensities'].astype(str).apply(lambda x: x[1:-1])
	df['chemical_formula'] = df['chemical_formula'].astype(str)

	# Clean up Intensity data for normalization
	to_delete = string.punctuation.replace('.', '')  # remove fullstop

	# Normalize intensities
	for row in df.itertuples():  # Normalise the intensity list by sum
		# this map may need to happen for unlabeled (Day 0) intensities
		temp = list(map(
			float, row.intensities.translate(
				str.maketrans('', '', to_delete)).split()
		))
		s = sum(temp)
		df.at[row.Index, 'intensities'] = [float(i) / s for i in temp]

	df = df.reset_index(drop=False)
	df = df[['index', 'chemical_formula', 'intensities', 'enrichment']]

	start = time.perf_counter()

	####
	#### Create Multiple Calculators Based on Chunk Size
	####

	chunkedData = []
	chunkSize = 1500
	numChunks = 0

	# Creates a temporary director to store intermediate files before combining to all for multiple instances running at the same time.
	# temp_dir = tempfile.mkdtemp()

	# Create a folder to output intermediate data to
	# inter_dir = (os.getcwd() + "/Intermediates")
	# if (not os.path.exists(inter_dir)):
	# 	os.mkdir(inter_dir)
	# inter_dir = (os.getcwd() + "/Intermediates/Intermediates_" + str(filename)[filename.rfind('/') + 1:-4])
	# print("Saving Intermediate files at: ", inter_dir)
	# if (os.path.exists(inter_dir)):
	# 	shutil.rmtree(inter_dir)
	# os.mkdir(inter_dir)

	# Create the Chunks
	for num in range(int(df.shape[0] / chunkSize)):
		chunkedData.append(
			df[(numChunks * chunkSize):((numChunks + 1) * chunkSize)])
		numChunks += 1
	if (numChunks < (df.shape[0] / chunkSize)):
		chunkedData.append(df[(numChunks * chunkSize):])
		numChunks += 1
	print("Size of each Chunk: " + str(chunkSize))
	print("Number of Rows to process: " + str(df.shape[0]))

	# Run the Chunks
	chunked_df = list()

	currentChunk = 1
	print("There are " + str(numChunks) + " chunk(s) to compute")
	for chunk in chunkedData:
		# Create the Calculators
		calculator = NValueCalculator(
			chunk, output_columns=output_columns)
		print('')
		print(f"***Currently Running Chunk: " + str(currentChunk) + f"/" + str(numChunks) + f"*** ", end="\r",
		      flush=True)

		calculator.run()
		calculator.df = calculator.df.sort_index()
		chunked_df.append(calculator.df)
		currentChunk += 1
		del calculator  # Free Memory

	combined_csv = pd.concat(chunked_df)  # combine all files in the list

	# Merge Intermediates with original file
	combined_csv = combined_csv.set_index('index')
	combined_csv = combined_csv.sort_index()
	df_full = df_full.merge(
		right=combined_csv[output_columns],
		how='outer',
		left_index=True,
		right_index=True
	)

	# Output final files and show time elapsed
	end = time.perf_counter()

	print('')
	print('Process took {} seconds to run'.format(end - start))

	return df_full

def main():
	root = Tk()
	root.withdraw()
	root.filename = filedialog.askopenfilename(
		initialdir=os.getcwd(), title="Select file", filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
	root.update()

	filename = root.filename

	results = data_chunker(pd.read_csv(filename))

	print("hello")

if __name__ == "__main__":
	main()