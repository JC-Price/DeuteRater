from utils.NValueCalculator import NValueCalculator
import time
import pandas as pd
import string

from tkinter import filedialog
from tkinter import Tk
import os

output_columns = ['empir_n', 'stddev']

def data_chunker(df_full):
	# Renames columns to proper names depending on if a lipid or protein is being passed in
	df = df_full.rename(columns={  # TODO: Verify Column Names
		'cf': 'chemical_formula',
		'abundances': 'intensities',
		'labeling': 'enrichment'})

	df['intensities'] = df['intensities'].astype(str).apply(lambda x: x[1:-1])
	df['chemical_formula'] = df['chemical_formula'].astype(str)

	# Normalize intensities
	for row in df.itertuples():  # Normalise the intensity list by sum
		temp = [float(x) for x in row.intensities.split(', ')]
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
	chunked_dfs = list()

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
		chunked_dfs.append(calculator.df)
		currentChunk += 1
		del calculator  # Free Memory

	combined_chunks = pd.concat(chunked_dfs)  # combine all files in the list

	# Merge Intermediates with original file
	combined_chunks = combined_chunks.set_index('index')
	combined_chunks = combined_chunks.sort_index()
	df_full = df_full.merge(
		right=combined_chunks[output_columns],
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