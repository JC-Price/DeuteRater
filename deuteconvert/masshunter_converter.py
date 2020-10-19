# This code was made to take in 3 files from PEAKs:
# feature, peptides, protein-peptides
import pandas as pd
import numpy as np
from functools import partial
import json
import re

from base_converter import BaseConverter

import peptide_utils as peputils
import settings


class MassHunterConverter(BaseConverter):
    # TODO: slots

    settings.load('../resources/converter_settings.yaml')
    # TODO: These constants need organized as well
    accession_parse_token = '|'
    PROTON_MASS = 1.007825

    # TODO: determine which columns are required

    # These are the most likely things to change in different versions
    correct_header_names = {
        'Name': 'Lipid Name',
        'Formula': 'cf',
        'Score': 'Score',
        'Mass': 'theoretical_mass',
        'm/z': 'Precursor m/z',

        # These numbers don't look right in the data. Check what they mean
        'RT': 'Precursor Retention Time (sec)',
        'Start': 'rt_start',
        'End': 'rt_end',
        'Width': 'rt_width',

        # change the name of rt_mean to P.R.T (sec)
        # 'theoretical_mass': 'Peptide Theoretical Mass',
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

        'neutromers_to_extract': 'neutromers_to_extract',
        'literature_n': 'literature_n'
    }

    correct_header_order = [
        # Columns to Pass Through
        'Lipid Name',  # -> 'Lipid Name'
        'Hits',
        'HMP',
        'LMP',
        'Notes',
        'MH Abundance',  # -> 'MH Abundance'
        'Mining Algorithm',
        'Area',
        'Base Peak',
        'Z Count',
        'MS/MS Count',
        'ID File Precursor',  # -> 'ID File Precursor'
        'rt_end',  # -> 'rt_end'
        'cf',  # -> 'cf'
        'Height',
        'Ions',
        'Polarity',
        'Lipid Unique Identifier',  # -> 'Lipid Unique Identifier'
        'Min Z',
        'Max Z',
        'Score',  # (Equivalent to Quality)
        'Precursor m/z',  # -> 'Precursor m/z'
        'Precursor Retention Time (sec)',  # -> 'Precursor Retention Time (sec)'
        'Saturated',
        'rt_start',  # -> 'rt_start'
        'rt_width',  # -> 'rt_width'
        'Diff (Tgt, mDa)',
        'Diff (Tgt, ppm)',
        'Score (Tgt)',
        'Flags (Tgt)',
        'Flag Severity (Tgt)',
        'Mass (Tgt)',
        'Sample Name',
        'Instrument Name',
        'Position',
        'Acq Method',
        'Identification Charge',
        # #Old Names
        # 'z',
        # 'ptm',
        # 'avg_ppm',
        # 'start_loc',
        # 'end_loc',
        # 'num_peptides',
        # 'num_unique',
        # 'accessions',  # -> LMP
        # 'species',
        # 'gene_name',
        # 'protein_existence',
        # 'sequence_version',
        'neutromers_to_extract',
        'literature_n'
    ]

    def __init__(self, precursor):
        self.precursor_path = precursor
        self._id_df = None

        super().__init__()

    @property
    def converted_data(self):
        if self._id_df.size == 0:
            raise RuntimeError('The converter has not run to completion')
        return self._id_df

    def load_files(self, parameterslist):
        pass

    def convert(self):
        self._id_df = self._load_precursor()

        self._set_n_peaks()
        self._finalize()

    def write(self, out):
        self._id_df.to_csv(out, sep='\t', index=False)

    # NOTE: add functions for additional formats here if we need them

    def _load_precursor(self):
        df = pd.read_csv(self.precursor_path, skiprows=2, index_col=False)
        keep_cols = [
            'Name',  # -> 'Lipid Name'
            'Hits',
            'HMP',
            'LMP',
            'Notes',
            'Abund',  # -> 'MH Abundance'
            'Mining Algorithm',
            'Area',
            'Base Peak',
            'Z Count',
            'MS/MS Count',
            'File',  # -> 'ID File Precursor'
            'End',  # -> 'rt_end'
            'Formula',  # -> 'cf'
            'Height',
            'Ions',
            'Polarity',
            'Label',  # -> 'Lipid Unique Identifier'
            'Min Z',
            'Max Z',
            'Mass',
            'Score',  # (Equivalent to Quality)
            'm/z',  # -> 'Precursor m/z'
            'RT',  # -> 'Precursor Retention Time (sec)'
            'Saturated',
            'Start',  # -> 'rt_start'
            'Width',  # -> 'rt_width'
            'Diff (Tgt, mDa)',
            'Diff (Tgt, ppm)',
            'Score (Tgt)',
            'Flags (Tgt)',
            'Flag Severity (Tgt)',
            'Mass (Tgt)',
            'Sample Name',
            'Instrument Name',
            'Position',
            'Acq Method',
        ]
        rename_cols = {
            'Name': 'Lipid Name',
            'Abund': 'MH Abundance',
            'File': 'ID File Precursor',
            'End': 'rt_end',
            'Formula': 'cf',
            'Label': 'Lipid Unique Identifier',
            'm/z': 'Precursor m/z',
            'RT': 'Precursor Retention Time (sec)',
            'Start': 'rt_start',
            'Width': 'rt_width',
            'Mass': 'theoretical_mass'
        }
        dtype_dict = {
            'rt_start': np.float,
            'rt_end': np.float,
            'Precursor Retention Time (sec)': np.float,
            'theoretical_mass': np.float,
            'Precursor m/z': np.float
        }

        df = df[keep_cols].rename(columns=rename_cols)
        df = df.astype(dtype=dtype_dict)
        df['cf'] = df['cf'].astype(str)
        df['cf'] = df['cf'].str.replace(' ', '')
        df = df.dropna(subset=['Area'])
        df['Identification Charge'] = (df['theoretical_mass'] / df['Precursor m/z']).round(0).astype(int)

        return df

    def _set_n_peaks(self):
        for mass_cutoff, n_peaks in settings.mass_cutoffs:
            self._id_df.loc[
                self._id_df['theoretical_mass'] > mass_cutoff,
                'neutromers_to_extract'
            ] = n_peaks

    def _finalize(self):
        self._id_df.sort_index(inplace=True)
        self._id_df['literature_n'] = np.nan
        self._id_df = self._id_df[MassHunterConverter.correct_header_order].rename(
            columns=MassHunterConverter.correct_header_names
        )

def main():
    converter = MassHunterConverter("/mnt/wwn-0x5000c500c4ea04a1-part2/CQ1/CQ1049/Test_Data/sample_M_RD-3-009_A2D0M1_MSMS_rep1Q1-27.csv")
    converter.convert()
    converter.write("/mnt/wwn-0x5000c500c4ea04a1-part2/CQ1/CQ1049/Test_Data/sample_M_RD-3-009_A2D0M1_MSMS_rep1Q1-27_IDFILE.csv")

if __name__ == "__main__":
    main()
