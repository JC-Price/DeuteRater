# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
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

# This code was made to take in 3 files from PEAKs:
# feature, peptides, protein-peptides
# TODO: Drop duplicates after sorting for best results

import pandas as pd
import numpy as np
from functools import partial
import json
import re
import os

from copy import copy

from deuteconvert.base_converter import BaseConverter
import deuteconvert.peptide_utils as peputils
import deuteconvert.settings as settings


#location = os.path.abspath(sys.executable)
location = os.path.dirname(os.path.abspath(__file__))

main_location = os.path.dirname(location)
json_path = os.path.join(main_location, "resources", "ptms.json")


#$adjusted from peaks 8.5 converter for use with Peaks X+
class PeaksXplus(BaseConverter):
    accession_parse_token = '|'
    PROTON_MASS = 1.007825

    #$change the headers to ensure consistent naming
    correct_header_names = {
        'peptide': 'Sequence',
        'rt_mean': 'Precursor Retention Time (sec)',
        'rt_start': 'rt_start',
        'rt_end': 'rt_end',
        'rt_width': 'rt_width',
        'mz': 'Precursor m/z',
        'theoretical_mass': 'Peptide Theoretical Mass',
        'z': 'Identification Charge',
        'first_accession': 'Protein ID',
        'accessions': 'Homologous Proteins',
        'ptm': 'ptm',
        'avg_ppm': 'avg_ppm',
        'start_loc': 'start_loc',
        'end_loc': 'end_loc',
        'num_peptides': 'num_peptides',
        'num_unique': 'num_unique',
        'protein': 'Protein Name',
        'species': 'species',
        'gene_name': 'gene_name',
        'protein_existence': 'protein_existence',
        'sequence_version': 'sequence_version',
        'cf': 'cf',
        'neutromers_to_extract': 'neutromers_to_extract',
        'literature_n': 'literature_n'
    }
    #$ ensure consistent order
    correct_header_order = [
        'peptide',
        'first_accession',
        'protein',
        'rt_mean',
        'rt_start',
        'rt_end',
        'rt_width',
        'mz',
        'theoretical_mass',
        'z',
        'ptm',
        'avg_ppm',
        'start_loc',
        'end_loc',
        'num_peptides',
        'num_unique',
        'accessions',
        'species',
        'gene_name',
        'protein_existence',
        'sequence_version',
        'cf',
        'neutromers_to_extract',
        'literature_n'
    ]

    def __init__(self, in_files,
                 settings_path):
        self.prot_path = in_files[0]
        self.protpep_path = in_files[1]
        self.feat_path = in_files[2]
        self._id_df = None

        super().__init__()
        settings.load(settings_path)


    #$get the data and do the analysis
    def convert(self):
        df_proteins = self._prepare_proteins()
        df_protein_peptides = self._prepare_protein_peptides()
        df_feature = self._prepare_feature()
        self._id_df = PeaksXplus._merge_data(
            df_feature, df_protein_peptides, df_proteins
        )
        #$stps to do in order.  yes proximity filter should be in their twice
        #$first to save time, and second to ensure no duplicates have popped up in the 
        #$expand_to_Charge_States function.
        funcs = [
            PeaksXplus._interpret_aa_sequences,
            PeaksXplus._set_n_peaks,
            PeaksXplus._proximity_filter,
            PeaksXplus._finalize,
            PeaksXplus._expand_to_charge_states,
            PeaksXplus._proximity_filter
        ]
        for fn in funcs:
            self._id_df = fn(self._id_df)

    def write(self, out):
        self._id_df.to_csv(out, sep=',', index=False)

    @staticmethod
    def _parse_accession(acc):
        return acc.split(PeaksXplus.accession_parse_token)[1]
    
    #$Peaks has a unique way of writing ptms, need to get a dictionary of valid symbols
    #$ and check for mutations
    @staticmethod
    def _reversed_ptm(peptide_str):
        with open(json_path) as ptms_json: #open(settings.ptms_path) as ptms_json:
            known_ptms = json.load(ptms_json)
        reversed_dict = {}
        for key in known_ptms:
            for mini_list in known_ptms[key]:
                reversed_dict[mini_list[1]] = key
        
        all_modifications = []
        for key in reversed_dict:
            
            if key in peptide_str:
                all_modifications.append(reversed_dict[key])
        #$need to add mutations. inserts, substitutions, and deletions
        #$since substitutions can happen for all can't just add, so just
        #$look for a few keywords
        for mutation in ["(ins", "(sub", "(del"]:
            if mutation in peptide_str:
                all_modifications.append('Mutation')
                break
        
        return "; ".join(all_modifications)

    #$need to actually replace the modifications with the correct single letter code    
    @staticmethod
    def _process_ptm(ptm_string, peptide_sequence):
        seq = peptide_sequence
        ptms = [ptm.strip() for ptm in ptm_string.split(';')]
        # TODO: dynamically load ptm file reference? should use gui or file?
        with open(json_path) as ptms_json: #open(settings.ptms_path) as ptms_json:
            known_ptms = json.load(ptms_json)
        for ptm in ptms:
            if ptm in known_ptms.keys():
                for mod_pair in known_ptms[ptm]:
                    seq = seq.replace(mod_pair[1], mod_pair[0])
            elif ptm == 'Mutation':
                seq = re.sub(r'\(sub [a-zA-Z0-9]\)', '', seq)
                seq = re.sub(r'\(del [a-zA-Z0-9]\)', '', seq)
                seq = re.sub(r'\(ins\)', '', seq)
        return seq
    
    #$actually merge the three separate inputs
    @staticmethod
    def _merge_data(df_feature, df_protein_peptides, df_proteins):
        data = df_feature.merge(
            how='inner',
            right=df_protein_peptides,
            left_on='peptide',
            right_on='peptide'
        ).sort_values('first_accession')
        data = data.merge(
            how='inner',
            right=df_proteins,
            left_on='first_accession',
            right_on='accession'
        )
        data = data.drop('accession', axis=1)
        data = data.sort_values(
            by=['peptide', 'num_unique', 'rt_mean'],
            ascending=[True, False, True]
        ).reset_index(drop=True)
        return data

    #$deal with the peptide_features file
    def _prepare_feature(self):  # change RT mean (in feature) to rt_mean
        df = pd.read_csv(self.feat_path)
        
        #$I have not seen evidence for this X+, it is likley a possibility
        #$however, I can't uset RT for this and checking for extras is not 
        #$reasonable since RT begin and end also have RT. For now we will
        #$comment it out until we find out if these are included
        # Get retention time column names
        #rt_names = [col for col in df.columns if 'RT mean' in col]
        #del rt_names[-1]

        keep_cols = [
            'DB Peptide',
            'Denovo Peptide',
            'RT Begin',
            'RT End',
            'RT',
            'm/z', 'z',
            'Accession'
        ]
        #keep_cols.extend(rt_names)
        rename_cols = {
            'RT Begin': 'rt_start',
            'RT End': 'rt_end',
            'RT': 'rt_mean',
            'm/z': 'mz',
            'z': 'z',
            'Accession': 'accessions'
        }
        dtype_dict = {
            'rt_start': np.float,
            'rt_end': np.float,
            'rt_mean': np.float,
            'mz': np.float,
            'z': np.int8
        }

        df = df[keep_cols].rename(columns=rename_cols)
        df['peptide'] = df['DB Peptide'].fillna(df['Denovo Peptide'])
        df = df.drop(['DB Peptide', 'Denovo Peptide'], axis=1)
        df = df.astype(dtype=dtype_dict)
        df[['rt_start', 'rt_end', 'rt_mean']] = \
            df[['rt_start', 'rt_end', 'rt_mean']] * 60
        df['rt_width'] = df['rt_end'].values - df['rt_start'].values
        df['accessions'] = df['accessions'].str.split(':').fillna('').map(
            partial(map, PeaksXplus._parse_accession)).map(list)
        #$this is a bit of a change in that features annoyingly does not have
        #$the ptms marked.  it does have the notes in the sequence though
        #$so we'll make the ptm column. we'll try and just run it backwards
        for row in df.itertuples():
            df.at[row.Index, 'ptm'] = PeaksXplus._reversed_ptm(row.peptide)
        
        ptm_mask = (df['ptm'] != '')
        for row in df[ptm_mask].itertuples():
            df.at[row.Index, 'peptide'] = \
                PeaksXplus._process_ptm(row.ptm, row.peptide)
        df['first_accession'] = df['accessions'].str[0]

#####################################################
        # replace dashes (NaN values) with 0
        #for col in rt_names:
        #    df[col] = df[col].replace('-', 0)
        #$replacement for above.  comment out if when finish
        df['rt_mean'] = df['rt_mean'].replace('-', 0)


        # filter based on starting time
        df = df[df['rt_mean'] >= settings.start_time]

        # count number of reps in retention time
        #$ when figure out how reps are represented uncomment the loop
        df['num_reps'] = 0
        #for column_name in rt_names:
        #    df.update(df.loc[df[column_name] != 0]['num_reps'] + 1)

        # fix the column types
        #df = df.astype(dict.fromkeys(rt_names, np.float))

        # calculate the standard error of the mean
        """
        get rid of this block comment when we figure out the the reps in 
        X+ and remove the line following this comment
        df['std_error'] = df[rt_names].sem(axis='columns').abs()
        """
        df['std_error'] = 0
        df = df.astype({'std_error': np.float})
        df = df.sort_values(
            by=['std_error', 'num_reps', 'rt_width', 'z'],
            ascending=[True, False, True, True]
        )
        df = df.drop_duplicates(subset='peptide', keep='first')
#####################################################

        df = df.sort_values(by=['peptide'])
        return df
    
    #$prepare the protein peptide file, mostly renaming, trimming and some parsing.
    def _prepare_protein_peptides(self):
        df = pd.read_csv(self.protpep_path)
        #$quality has been moved around. of these files only in the features.csv
        #$which is constant currently.  conveniently we're not using that currently
        #$uncomment if necessary
        keep_cols = [
            'Peptide',
            'ppm',
            'Start',
            'End',
            'PTM'
        ]
        rename_cols = {
            'Peptide': 'peptide',
            'ppm': 'avg_ppm',
            'Start': 'start_loc',
            'End': 'end_loc',
            'PTM': 'ptm'
        }
        df = df[keep_cols].rename(columns=rename_cols)
        df = df.sort_values(by=['peptide'])
        df = df.drop_duplicates(subset='peptide', keep='first')
        df['ptm'] = df['ptm'].fillna('')
        ptm_mask = (df['ptm'] != '')
        for row in df[ptm_mask].itertuples():
            df.at[row.Index, 'peptide'] = \
                PeaksXplus._process_ptm(row.ptm, row.peptide)
        df['peptide'] = df['peptide'].str.split('.').str[1]
        df = df.drop('ptm', axis=1)
        return df
    #$prepare the protein file, mostly renaming, trimming and some parsing.
    def _prepare_proteins(self):
        df = pd.read_csv(self.prot_path)
        keep_cols = ['Accession', '#Peptides', '#Unique', 'Description']
        rename_cols = {
            'Accession': 'accession',
            '#Peptides': 'num_peptides',
            '#Unique': 'num_unique',
            'Description': 'protein'
        }
        df = df[keep_cols].rename(columns=rename_cols)
        df['accession'] = df['accession'].map(PeaksXplus._parse_accession)
        description_regex = r'(.*?)OS=(.*?)(?:GN=(.*?))?PE=(.*?)SV=(.*)'
        description_strings = df['protein'].str.extract(
            description_regex, expand=True
        )
        df['protein'] = description_strings[0]
        df['species'] = description_strings[1]
        df['gene_name'] = description_strings[2]
        df['protein_existence'] = description_strings[3]
        df['sequence_version'] = description_strings[4]
        return df

    #$turn the amino acids sequences into something we can calculate off of
    @staticmethod
    def _interpret_aa_sequences(df):
        df.assign(
            cf='',
            theoretical_mass=np.nan,
            literature_n=np.nan
        )
        aa_comp_df = pd.read_csv(settings.aa_elem_comp_path, sep='\t')
        aa_comp_df.set_index('amino_acid', inplace=True)

        aa_label_df = pd.read_csv(settings.aa_label_path, sep='\t')
        aa_label_df.set_index('study_type', inplace=True)
        # TODO: Find out where to store settings, then decide which 'studytype'
        #       to use as default.
        aa_labeling_dict = aa_label_df.loc[settings.study_type, ].to_dict()

        elem_df = pd.read_csv(settings.elems_path, sep='\t')
        # This is necessary if we have all of the different isotopes in the tsv
        element_index_mask = [
            0,  # Hydrogen
            10,  # Carbon-12
            13,  # Nitrogen-14
            15,  # Oxygen-16,
            30, # Phosphorous
            31  # Sulfer-32
        ]
        elem_df = elem_df.iloc[element_index_mask]
        elem_df.set_index('isotope_letter', inplace=True)

        for row in df.itertuples():
            i = row.Index
            aa_counts = {}
            for aa in row.peptide:
                if aa not in aa_counts.keys():
                    aa_counts[aa] = 0
                aa_counts[aa] += 1
            elem_dict = peputils.calc_cf(aa_counts, aa_comp_df)
            theoretical_mass = peputils.calc_theory_mass(elem_dict, elem_df)
            literature_n = peputils.calc_add_n(aa_counts, aa_labeling_dict)

            df.at[i, 'cf'] = ''.join(
                k + str(v) for k, v in elem_dict.items() if v > 0
            )
            df.at[i, 'theoretical_mass'] = theoretical_mass
            df.at[i, 'literature_n'] = literature_n
        return df

    #$determine the number of isotope peaks to extract
    @staticmethod
    def _set_n_peaks(df):
        for mass_cutoff, n_peaks in settings.mass_cutoffs:
            df.loc[
                df['theoretical_mass'] > mass_cutoff,
                'neutromers_to_extract'
            ] = n_peaks
        return df

    #$ppm calculation
    @staticmethod
    def _ppm_calculator(target, actual):
        ppm = (target-actual)/target * 1000000
        return abs(ppm)
    #$since we are extracting based on mz and retention time we need
    #$to ensure they are unique.  this removes non-unique retention times and m/z pairs
    #$ based on settings the user provides
    @staticmethod
    def _proximity_filter(df):
        try:
            df.sort_values(by='mz', inplace=True)
            mz_index = list(df.columns).index("mz") 
            rt_index = list(df.columns).index("rt_mean")
        except KeyError:
            df.sort_values(by='Precursor m/z', inplace=True)
            mz_index = list(df.columns).index('Precursor m/z') 
            rt_index = list(df.columns).index("Precursor Retention Time (sec)")
        list_of_lists = df.values.tolist()
        too_close = []
        for i in range(len(list_of_lists)):
            for j in range(i+1, len(list_of_lists)):
                current_ppm = PeaksXplus._ppm_calculator(list_of_lists[i][mz_index], list_of_lists[j][mz_index])
                if current_ppm > settings.mz_proximity_tolerance: 
                    break
                if abs(list_of_lists[i][rt_index] - list_of_lists[j][rt_index]) < settings.rt_proximity_tolerance:
                    too_close.extend([i,j])
        too_close = list(set(too_close))
        return df.drop(df.index[too_close])

    #$expand to multiple charge states that may be present but were not identified
    @staticmethod
    def _expand_to_charge_states(df):
        all_columns = list(df.columns)
        precursor_mz_index = list(df.columns).index('Precursor m/z') 
        id_charge_index = list(df.columns).index('Identification Charge')
        list_of_lists = df.values.tolist()
        all_data = []
        for sub_list in list_of_lists:
            premass = sub_list[precursor_mz_index] * sub_list[id_charge_index] - \
                sub_list[id_charge_index] * PeaksXplus.PROTON_MASS
            for z in range(settings.min_charge_state, settings.max_charge_state + 1):
                quick_list =copy(sub_list)
                quick_list[precursor_mz_index] = (premass + float(z) * PeaksXplus.PROTON_MASS) / float(z)
                quick_list[id_charge_index] = z
                all_data.append(quick_list)
        return pd.DataFrame(all_data, columns = all_columns)

    #$need to sort, rename columns and mark too short peptides.
    @staticmethod
    def _finalize(df):
        df = df[df['peptide'].str.len() >= 6]
        df = df.sort_index()
        df = df[PeaksXplus.correct_header_order].rename(
            columns=PeaksXplus.correct_header_names
        )
        return df
