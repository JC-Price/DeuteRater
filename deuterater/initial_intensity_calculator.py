# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, Benjamin Driggs,
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
import traceback

"""
as with combine_extracted_files, this provides some basic filtering and 
calculations that are relatively simple on their own, but need to be done and it's
easier here than in the main calculation function

The main filters here are n value and sequence length must be sufficiently large
to be confident we are going to actually see good signal to noise. Finally we must
use the emass algorithm to calculate the unenriched relative isotopic envelope

may merge into combine_extracted_files.py
"""

from tqdm import tqdm
import pandas as pd
import multiprocessing as mp
import concurrent.futures as cf
import numpy as np

from pathlib import Path
from functools import partial

from utils.emass import emass
import deuterater.settings as settings
import utils.old_NValueCalculator as nvct
import resources.peptide_utils as peptide_utils

max_isos = 5  # constant based on the n_isos based on the mass (done in the extractor)
p0_guess = 1  # seems to work for most fits. if it causes problems we can adjust

# used in the _calculate_literature_n() function brought over from non-human version of DeuteRater
literature_n_name = "literature_n"
sequence_column_name = "Sequence"


# as with all the calculation steps this is a class for consistent calls in the main
class theoretical_enrichment_calculator(object):
    def __init__(self, prepared_data_path, out_path, settings_path, biomolecule_type):
        settings.load(settings_path)
        self.settings_path = settings_path

        self.prepared_data_path = Path(prepared_data_path)
        self.out_path = out_path
        self.biomolecule_type = biomolecule_type

        if self.biomolecule_type == "Peptide":
            aa_label_df = pd.read_csv(settings.aa_labeling_sites_path, sep='\t')
            aa_label_df.set_index('study_type', inplace=True)
            self.aa_labeling_dict = aa_label_df.loc[settings.label_key,].to_dict()
            # settings.use_empir_n_value = False

        if self.prepared_data_path.suffix == '.tsv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep='\t'
            )
        elif self.prepared_data_path.suffix == '.csv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep=','
            )
        # if multiprocessing need to set that up. more than 60 cores causes problems for windows
        if settings.recognize_available_cores is True:
            # BD: Issue with mp.cpu_count() finding too many cores available
            self._n_processors = round(mp.cpu_count() * 0.80)
        else:
            self._n_processors = settings.n_processors
        if self._n_processors > 60:
            self.n_processors = 60
        self.model = None

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    # Run the analysis. This function doesn't have any calculation itself (other than merging, transposing,
    # and setting index of the results) it prepares a function for multiprocessing and then begins the multiprocessing
    def prepare(self):
        if self.biomolecule_type == "Peptide":
            unique_molecules_df = self.data_df.drop_duplicates(subset=["Sequence"])
        else:
            unique_molecules_df = self.data_df.drop_duplicates()

        new_columns = _make_new_columns(self.biomolecule_type)
        func = partial(theoretical_enrichment_calculator._individual_process,
                       new_columns=new_columns,
                       minimum_n_value=settings.min_allowed_n_values,
                       minimum_sequence_length=settings.min_aa_sequence_length,
                       biomolecule_type=self.biomolecule_type,
                       use_empir_n_value=settings.use_empir_n_value)
        try:
            df_split = np.array_split(unique_molecules_df, len(unique_molecules_df))
        except Exception as e:
            print("ERROR: No valid data in combined_extracted_files output file")
            print(e)
            traceback.print_tb(e.__traceback__)
            raise

        # TODO: add a debug version of this process. - Ben D

        if settings.debug_level == 0:
            with cf.ProcessPoolExecutor(max_workers=self._n_processors) as executor:
                final_df = pd.concat(
                    tqdm(executor.map(func, df_split), total=len(df_split), desc="Calculating Initial Intensities: "),
                    axis=1).drop_duplicates(keep='first')
        elif settings.debug_level >= 1:
            dfs = []
            for dframe in tqdm(df_split, desc="Calculating Initial Intensities: "):
                dfs.append(func(dframe))
            final_df = pd.concat(dfs, axis=1).drop_duplicates(keep='first')

        final_df = final_df.T
        if self.biomolecule_type == "Peptide":
            final_df = final_df.set_index("Sequence")
            self.model = pd.merge(self.data_df, final_df, left_on="Sequence", right_index=True).drop_duplicates()
        else:
            final_df = final_df.set_index("Adduct_cf")
            self.model = pd.merge(self.data_df, final_df, left_on="Adduct_cf", right_index=True).drop_duplicates()

    # actually runs the relevant calculation.
    @staticmethod
    def _individual_process(df, new_columns, minimum_n_value, minimum_sequence_length, biomolecule_type,
                            use_empir_n_value):
        variable_list = []
        for row in df.itertuples():
            output_series = pd.Series(index=new_columns, dtype="object")

            if biomolecule_type == "Peptide":
                output_series["Sequence"] = row.Sequence
            else:
                output_series["Adduct_cf"] = row.Adduct_cf

            # check if there is a valid n_value, otherwise we'll drop the row
            if use_empir_n_value:
                if row.n_value == "no valid time points" or row.n_value == "error occurred":
                    variable_list.append(_error_message_results("Error: see n_value column", output_series))
                    continue

            try:
                # drop rows we will not use before doing any compiles calculations
                if biomolecule_type == "Peptide" and len(row.Sequence) < minimum_sequence_length:
                    variable_list.append(_error_message_results(
                        f"Sequence is less than {minimum_sequence_length} amino acids", output_series))
                    continue
                if biomolecule_type == "Lipid" and float(row.n_value) < minimum_n_value:
                    variable_list.append(_error_message_results(
                        f"less than {minimum_n_value} labeling sites",
                        output_series))
                    continue
            except:
                # n-value column already has error message, so we just continue with the calculations
                continue

            intensity_values = _fit_emass(row.cf, row.n_isos)

            output_series["Theoretical Unlabeled Normalized Abundances"] = ", ".join(intensity_values)
            variable_list.append(output_series)
        return pd.concat(variable_list, axis=1)


# if an error happens it is most efficient to have a function
def _error_message_results(error_message, output_series):
    # don't need to know which names are which or how many columns there are, 
    # just need python to fill all non-Sequence columns
    # position 0 is the sequence name which we don't wish to overwrite
    for index_name in output_series.index[1:]:
        output_series[index_name] = error_message
    return output_series


# calculate unlabeled intensity, if we need to return  m/z values or
# adjust for different n_values, do it here or in emass itself.
def _fit_emass(cf, n_isos):
    intensity_values = emass(cf, n_isos)
    return [str(i) for i in intensity_values]


# this creates the header for the variables. is a function in case we need 
# to add various columns  (if we want to graph emass output or something)
def _make_new_columns(biomolecule_type):
    if biomolecule_type == "Peptide":
        new_columns = ["Sequence"]
    else:
        # TODO: What are these columns for and do we need them for lipids? - Ben D
        new_columns = ["Adduct_cf"]
    new_columns.extend(["Theoretical Unlabeled Normalized Abundances"])
    return new_columns
