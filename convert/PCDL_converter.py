# -*- coding: utf-8 -*-
"""

"""
# This code was made to take in 3 files from PEAKs:
import pandas as pd
import numpy as np
from copy import deepcopy

#from functools import partial
#import json
#import re

try:
    from .base_converter import BaseConverter
    from convert import settings
except:
    from base_converter import BaseConverter
    import settings

def tk_get_single_file(extension='*', prompt="Select file"):
    from tkinter import filedialog
    from tkinter import Tk
    import os
    if extension[0] != '*':
        extension = '*' + extension
    root = Tk()
    root.withdraw()
    if (extension == '*'):
        root.filename = filedialog.askopenfilename(
            initialdir=os.getcwd(), title=prompt)
    else:
        root.filename = filedialog.askopenfilename(
            initialdir=os.getcwd(), title=prompt,
            filetypes=((extension[2:] + " files", extension)))
    root.update()

    filename = root.filename
    return filename  # A string representing the file path

ELECTRON_MASS = 0.000549

atomic_mass =    {'H': 1.0078246,
                  'C': 12.0000000,
                  'N': 14.0030732,
                  'O': 15.9949141,
                  'P': 30.973762,
                  'S': 31.972070,
                  'F': 18.99840322,
                  'D': 2.0141021,
                  'Cl': 34.9688527,
                  'Br': 78.9183376,
                  'I': 126.9044719,
                  'Si': 27.9769265,
                  'Na': 22.98976928}

class PCDL_Converter(BaseConverter):
    # TODO: slots

    # settings.load('../resources/guide_settings.yaml')
    # TODO: These constants need organized as well
    accession_parse_token = '|'
    PROTON_MASS = 1.007825

    # TODO: determine which columns are required

    # These are the most likely things to change in different versions
    correct_header_names = {
         'RT': 'Precursor Retention Time (sec)',
         "Mass": 'Precursor m/z',
         "": 'Identification Charge',
         "": "Lipid Unique Identifier",
         "": "LMP",
         "Name ": "Lipid Name",
         "": "HMP",
         "Formula": "cf"
    }

    PCDL_correct_header_names = {
        "Lipid Name": "Name",
        "cf": "Formula",
        "Precursor m/z": "Mass",
        "Precursor Retention Time (min)": "Retention Time"
    }
    
    PCDL_correct_header_order = [
        "Name", "Formula", "Mass", "Retention Time", "Retention Index", "Cation", "Anion", "CAS", "ChemSpider",
        "PubChem", "Synonyms", "IUPAC", "NumSpectra", "CCS Count"
    ]

    PCDL_development_header_order = [
        "Name", "Formula", "Mass", "Retention Time", "Retention Index", "Cation", "Anion", "CAS", "ChemSpider",
        "PubChem", "Synonyms", "IUPAC", "NumSpectra", "CCS Count",
        'Matched_Results_Analysis', 'Matched_Details_Replicates_Used'
    ]

    correct_header_order = [
        'Lipid Name',
        'Lipid Unique Identifier',
        'Precursor m/z',
        'Precursor Retention Time (sec)',
        "Identification Charge",
        'LMP',
        'HMP',
        'cf',
        'neutromers_to_extract',
        'literature_n'
    ]

    development_header_order = [
        'Lipid Name',
        'Lipid Unique Identifier',
        'Precursor m/z',
        'Precursor Retention Time (sec)',
        "Identification Charge",
        'LMP',
        'HMP',
        'cf',
        'neutromers_to_extract',
        'literature_n',
        'Matched_Results_Analysis',
        'Matched_Details_Replicates_Used'
    ]

    def __init__(self, precursor, settings_path="../resources/id_settings.yaml"):
        self.precursor_path = precursor[0]
        self._id_df = None

        super().__init__()
        settings.load(settings_path)

    @property
    def converted_data(self):
        if self._id_df.size == 0:
            raise RuntimeError('The converter has not run to completion')
        return self._id_df

    def load_files(self, parameterslist):
        pass

    def convert(self):
        self._id_df = self._load_precursor()

        # if "Fragment" in self._id_df:
        #     self._id_df["Lipid Name"] = self._id_df["Lipid Name"] + self._id_df.Fragment

        self._set_n_peaks()
        # if settings.remove_duplicates:
        #     self._drop_duplicate_molecules()
        # self._id_df = self._proximity_filter(self._id_df)
        self._finalize()
        if settings.expand_to_adducts:
            self._id_df = self._expand_to_adducts(self._id_df)
        else:
            self._id_df = self._correct_adduct_cf(self._id_df)
        self._expand_to_charge_states()
        # self._id_df = self._proximity_filter(self._id_df)

    def write(self, out):
        self._id_df.to_csv(out, index=False)

    # NOTE: add functions for additional formats here if we need them

    def _drop_duplicate_molecules(self):
        df = self._id_df
        if "Fragment" in self._id_df.columns:
            self._id_df = df.loc[~df.duplicated(subset=["Lipid Name", "Fragment"], keep=False)]
        else:
            self._id_df = df.loc[~df.duplicated(subset="Lipid Name", keep=False)]
        print("Unique Lipids b4 duplicate removal = ", np.unique(df["Lipid Name"]).size)
        print("Unique Lipids after duplicate removal = ", np.unique(self._id_df["Lipid Name"]).size)

    def _load_precursor(self):
        df = pd.read_csv(self.precursor_path)
        rename_cols = {
            'Name': 'Lipid Name',
            'Formula': 'cf',
            'Mass': 'Precursor m/z',
            'Retention Time': 'Precursor Retention Time (sec)',
        }
        if "Ion Species" in df.columns:
            rename_cols["Ion Species"] = "Adduct"

        keep_cols = ['Name', 'Formula', 'Mass', 'Retention Time']

        dtype_dict = {
            'Precursor Retention Time (sec)': np.float64,
            'Precursor m/z': np.float64
        }

        # df = df[keep_cols].rename(columns=rename_cols)
        df = df.rename(columns=rename_cols)
        df = df.astype(dtype=dtype_dict)
        df['cf'] = df['cf'].astype(str)
        df['Identification Charge'] = 1
        df['Precursor Retention Time (min)'] = df["Precursor Retention Time (sec)"]
        df['Precursor Retention Time (sec)'] = df["Precursor Retention Time (sec)"] * 60

        return df

    @staticmethod
    def _ppm_calculator(target, actual):
        ppm = (target - actual) / target * 1000000
        return abs(ppm)

    @staticmethod
    def _expand_to_adducts(df):
        expanded_df = pd.DataFrame(columns=df.columns)
        for row in df.itertuples():
            try:
                temp_row = df.loc[row[0]]
            except:
                print("Oops, what happened")
                
            formmap = PCDL_Converter._parse_cf(row.cf)

            default_mz = temp_row['Precursor m/z']
            
            if settings.H_adduct:
                hydrogen_adduct = deepcopy(formmap)
                hydrogen_adduct['H'] += 1
                hydrogen_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in hydrogen_adduct.items())
                hydrogen_mz = default_mz + atomic_mass["H"]
                temp_row['Adduct'] = 'M+H'
                temp_row['Adduct_cf'] = hydrogen_cf
                temp_row['Precursor m/z'] = hydrogen_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.Na_adduct:
                sodium_adduct = deepcopy(formmap)
                sodium_adduct['Na'] += 1
                sodium_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in sodium_adduct.items())
                sodium_mz = default_mz + atomic_mass["Na"]
                temp_row['Adduct'] = 'M+Na'
                temp_row['Adduct_cf'] = sodium_cf
                temp_row['Precursor m/z'] = sodium_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.H_H2O_adduct:
                if "O" in formmap:
                    hydrogen_minus_water = deepcopy(formmap)
                    hydrogen_minus_water['H'] -= 1
                    hydrogen_minus_water['O'] -= 1
                    hydrogen_minus_water_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in hydrogen_minus_water.items())
                    hydrogen_minus_water_mz = default_mz - atomic_mass["H"] - atomic_mass["O"]
                    temp_row['Adduct'] = 'M+H-[H20]'
                    temp_row['Adduct_cf'] = hydrogen_minus_water_cf
                    temp_row['Precursor m/z'] = hydrogen_minus_water_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.Na_H2O_adduct:
                if "O" in formmap:
                    sodium_minus_water = deepcopy(formmap)
                    sodium_minus_water['Na'] += 1
                    sodium_minus_water['H'] -= 2
                    sodium_minus_water['O'] -= 1
                    sodium_minus_water_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in sodium_minus_water.items())
                    sodium_minus_water_mz = default_mz + atomic_mass["Na"] - atomic_mass["H"] * 2 - atomic_mass["O"]
                    temp_row['Adduct'] = 'M+Na-[H2O]'
                    temp_row['Adduct_cf'] = sodium_minus_water_cf
                    temp_row['Precursor m/z'] = sodium_minus_water_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.NH4_adduct:
                ammonium_adduct = deepcopy(formmap)
                ammonium_adduct['N'] += 1
                ammonium_adduct['H'] += 4
                ammonium_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in ammonium_adduct.items())
                ammonium_adduct_mz = default_mz + atomic_mass["N"] + atomic_mass["H"] * 4
                temp_row['Adduct'] = 'M+NH4'
                temp_row['Adduct_cf'] = ammonium_cf
                temp_row['Precursor m/z'] = ammonium_adduct_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)
            
            if settings.NH4_H2O_adduct:
                if "O" in formmap:
                    ammonium_minus_water = deepcopy(formmap)
                    ammonium_minus_water['N'] += 1
                    ammonium_minus_water['H'] += 2
                    ammonium_minus_water['O'] -= 1
                    ammonium_minus_water_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in ammonium_minus_water.items())
                    ammonium_minus_water_mz = default_mz + atomic_mass["N"] + atomic_mass["H"] - atomic_mass["O"]
                    temp_row['Adduct'] = 'M+NH4-[H2O]'
                    temp_row['Adduct_cf'] = ammonium_minus_water_cf
                    temp_row['Precursor m/z'] = ammonium_minus_water_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)
        
            if settings.formate_adduct:
                # Formate:
                formate_adduct = deepcopy(formmap)
                formate_adduct['C'] += 1
                formate_adduct['H'] += 1
                formate_adduct['O'] += 2
                formate_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in formate_adduct.items())
                formate_mz = default_mz + atomic_mass["H"] + atomic_mass["C"] + atomic_mass["O"] * 2 + ELECTRON_MASS
                temp_row['Adduct'] = '+COOH-'
                temp_row['Adduct_cf'] = formate_cf
                temp_row['Precursor m/z'] = formate_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)
            
            if settings.acetate_adduct:
                # Acetate:
                acetate_adduct = deepcopy(formmap)
                acetate_adduct['C'] += 2
                acetate_adduct['H'] += 3
                acetate_adduct['O'] += 2
                acetate_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in acetate_adduct.items())
                acetate_mz = default_mz + atomic_mass["H"] + atomic_mass["C"] + atomic_mass["O"] * 2 + ELECTRON_MASS
                temp_row['Adduct'] = '+C2H3O2-'
                temp_row['Adduct_cf'] = acetate_cf
                temp_row['Precursor m/z'] = acetate_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.e_adduct:
                # electron:
                electron_adduct = deepcopy(formmap)
                electron_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in electron_adduct.items())
                electron_mz = default_mz + ELECTRON_MASS
                temp_row['Adduct'] = '+e-'
                temp_row['Adduct_cf'] = electron_cf
                temp_row['Precursor m/z'] = electron_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.H_loss:
                # Hydrogen Loss:
                hydrogen_abduct = deepcopy(formmap)
                hydrogen_abduct['H'] -= 1
                hydrogen_abduct_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in hydrogen_abduct.items())
                hydrogen_abduct_mz = default_mz - atomic_mass["H"] + ELECTRON_MASS
                temp_row['Adduct'] = '-H'
                temp_row['Adduct_cf'] = hydrogen_abduct_cf
                temp_row['Precursor m/z'] = hydrogen_abduct_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.H2O_loss:
                if "O" in formmap:
                    # Water Loss:
                    water_abduct = deepcopy(formmap)
                    water_abduct['H'] -= 1
                    water_abduct_cf = ''.join(
                        '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in water_abduct.items())
                    water_abduct_mz = default_mz - atomic_mass["H"] * 2 - atomic_mass["O"] + ELECTRON_MASS
                    temp_row['Adduct'] = '-H2O'
                    temp_row['Adduct_cf'] = water_abduct_cf
                    temp_row['Precursor m/z'] = water_abduct_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)

        return expanded_df
    
    @staticmethod
    def _correct_adduct_cf(df):
        for idx, row in df.iterrows():
            formmap = PCDL_Converter._parse_cf(row.cf)
            if row.Adduct == "[M+H]+":
                hydrogen_adduct = deepcopy(formmap)
                hydrogen_adduct['H'] += 1
                hydrogen_cf = ''.join('%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in hydrogen_adduct.items())
                df.at[idx, 'Adduct'] = 'M+H'
                df.at[idx, 'Adduct_cf'] = hydrogen_cf
            elif row.Adduct == "[M+H4N]+":
                ammonium_adduct = deepcopy(formmap)
                ammonium_adduct['N'] += 1
                ammonium_adduct['H'] += 4
                ammonium_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in ammonium_adduct.items())
                df.at[idx, 'Adduct'] = 'M+NH4'
                df.at[idx, 'Adduct_cf'] = ammonium_cf
            elif row.Adduct == "[M-H]-":
                hydrogen_abduct = deepcopy(formmap)
                hydrogen_abduct['H'] -= 1
                hydrogen_abduct_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in hydrogen_abduct.items())
                df.at[idx, 'Adduct'] = '-H'
                df.at[idx, 'Adduct_cf'] = hydrogen_abduct_cf
            elif row.Adduct == "[M+CHO2]-":  # Formate:
                formate_adduct = deepcopy(formmap)
                formate_adduct['C'] += 1
                formate_adduct['H'] += 1
                formate_adduct['O'] += 2
                formate_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in formate_adduct.items())
                df.at[idx, 'Adduct'] = '+COOH-'
                df.at[idx, 'Adduct_cf'] = formate_cf
            elif row.Adduct == "[M+C2H3O2]-":  # Acetate:
                acetate_adduct = deepcopy(formmap)
                acetate_adduct['C'] += 2
                acetate_adduct['H'] += 3
                acetate_adduct['O'] += 2
                acetate_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % (k) if v == 1 else '' for k, v in acetate_adduct.items())
                df.at[idx, 'Adduct'] = '+C2H3O2-'
                df.at[idx, 'Adduct_cf'] = acetate_cf
            else:
                raise Exception(f"Unrecognized Adduct {row.Adduct}")
        return df
    
    @staticmethod
    def _parse_cf(cf):
        import re
        input_string = cf
        formmap = {}
        elements_separated = re.findall('[A-Z][^A-Z]*', input_string)
        for e in elements_separated:
            element = e.rstrip('0123456789')
            number = e[len(element):]
            if number == "": number = 1  # $ need the one to be there if element is alone
            formmap[element] = int(number)
    
        # add extra items if they don't exist
        if 'N' not in formmap:
            formmap['N'] = 0
        if 'Na' not in formmap:
            formmap['Na'] = 0
        if 'H' not in formmap:
            formmap['H'] = 0
        if 'O' not in formmap:
            formmap['O'] = 0
            
        return formmap

    @staticmethod
    def _proximity_filter(df):
        df.sort_values(by='Precursor m/z', inplace=True)
        mz_index = list(df.columns).index("Precursor m/z")
        rt_index = list(df.columns).index("Precursor Retention Time (sec)")
        list_of_lists = df.values.tolist()
        too_close = []
        for i in range(len(list_of_lists)):
            for j in range(i + 1, len(list_of_lists)):
                current_ppm = PCDL_Converter._ppm_calculator(list_of_lists[i][mz_index], list_of_lists[j][mz_index])
                if current_ppm > settings.mz_proximity_tolerance:
                    break
                # Multiply value in settings file to get seconds
                if abs(list_of_lists[i][rt_index] - list_of_lists[j][rt_index]) < settings.rt_proximity_tolerance * 60:
                    too_close.extend([i, j])
        too_close = list(set(too_close))
        return df.drop(df.index[too_close])

    def _set_n_peaks(self):
        # self._id_df['neutromers_to_extract'] = 3
        for mass_cutoff, n_peaks in settings.mass_cutoffs:
            self._id_df.loc[
                self._id_df['Precursor m/z'] > mass_cutoff,
                'neutromers_to_extract'
            ] = n_peaks

    def _expand_to_charge_states(self):
        from copy import copy
        df = self._id_df
        all_columns = list(df.columns)
        precursor_mz_index = list(df.columns).index('Precursor m/z')
        id_charge_index = list(df.columns).index('Identification Charge')
        list_of_lists = df.values.tolist()
        all_data = []
        for sub_list in list_of_lists:
            premass = sub_list[precursor_mz_index] * sub_list[id_charge_index] - \
                      sub_list[id_charge_index] * PCDL_Converter.PROTON_MASS
            for z in range(settings.min_charge_state, settings.max_charge_state + 1):
                quick_list = copy(sub_list)
                quick_list[precursor_mz_index] = (premass + float(z) * PCDL_Converter.PROTON_MASS) / float(z)
                quick_list[id_charge_index] = z
                all_data.append(quick_list)
        self._id_df = pd.DataFrame(all_data, columns=all_columns)

    def _finalize(self):
        header_order = PCDL_Converter.correct_header_order
        if "Fragment" in self._id_df.columns:
            PCDL_Converter.correct_header_order.insert(2, "Fragment")
        if settings.development:
            header_order = PCDL_Converter.development_header_order
            if "Extracted Comparison Results" in self._id_df.columns:
                header_order.append("Extracted Comparison Results")
            if "Concentration Dependence Results" in self._id_df.columns:
                header_order.append("Concentration Dependence Results")

        self._id_df.sort_index(inplace=True)
        self._id_df['literature_n'] = np.nan
        self._id_df['HMP'] = np.nan
        self._id_df['LMP'] = np.nan
        self._id_df['Lipid Unique Identifier'] = self._id_df['Lipid Name'] + "_" + self._id_df['Precursor Retention Time (sec)'].apply(lambda x: str("{:.3f}".format(x/60)))
        self._id_df = self._id_df[header_order].rename(
            columns=PCDL_Converter.correct_header_names
        )

def main():
    file = [tk_get_single_file()]
    converter = PCDL_Converter(file)
    converter.convert()
    print("done")

    converter.write(file[0][:-4] + "_ID_File.csv")

if __name__ == "__main__":
    main()
