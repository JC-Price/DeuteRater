# This code was made to take in 3 files from PEAKs:
import pandas as pd
import numpy as np
from copy import deepcopy

# from functools import partial
# import json
import re

try:
    from .base_converter import BaseConverter
    from . import settings
except:
    from base_converter import BaseConverter
    import settings


def tk_get_single_file(extension='*', prompt="Select file"):
    from tkinter import filedialog
    from tkinter import Tk
    import fnmatch
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

atomic_mass = {'H': 1.0078246,
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


class L4D_Analyzer_Converter(BaseConverter):
    # TODO: slots

    # TODO: These constants need organized as well
    accession_parse_token = '|'
    PROTON_MASS = 1.007825

    # TODO: determine which columns are required

    # These are the most likely things to change in different versions
    correct_header_names = {
        'CF': 'cf',
        'mz': 'Precursor m/z',
        'rt': 'Precursor Retention Time (sec)',
        'DT': 'Precursor Drift Time (ms)',
        'ccs': "CCS",
        'lipidMolecularSpecies': "Possible Annotations",
        'Lipid_Name': "Lipid Name",
    }

    correct_header_order = [
        "Lipid Name",
        "Lipid Unique Identifier",
        "Possible Annotations",
        'cf',
        'Adduct',
        'Adduct_cf',
        'Precursor m/z',
        'Precursor Retention Time (sec)',
        'Precursor Drift Time (ms)',
        'CCS',
        'Identification Charge',
        'neutromers_to_extract',
        "LMP",
        "HMP",
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

        self._set_n_peaks()
        self._id_df.dropna(subset=["Lipid Name"], inplace=True)
        self._remove_adduct_mass()
        self._id_df = self._expand_to_adducts(self._id_df)
        self._expand_to_charge_states()
        self._finalize()

    # remove unnamed lipids
    # self._id_df['dt'] = self._id_df["Precursor Drift Time (ms)"]
    # self._id_df['rt'] = self._id_df["Precursor Retention Time (sec)"] / 60
    # self._id_df['z'] = self._id_df["Identification Charge"]
    # self._id_df['mz'] = self._id_df["Precursor m/z"]
    # self._id_df['n_isos'] = self._id_df["neutromers_to_extract"]

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
        dtype_dict = {
            'Precursor Retention Time (sec)': np.float64,
            'Precursor m/z': np.float64,
            'Precursor Drift Time (ms)': np.float64,
            'cf': str,
        }

        df = pd.read_csv(self.precursor_path, sep='\t')

        # df = df[keep_cols].rename(columns=rename_cols)
        df = df.rename(columns=self.correct_header_names)
        df = df.astype(dtype=dtype_dict)
        df['Identification Charge'] = 1
        df['Precursor Retention Time (min)'] = df["Precursor Retention Time (sec)"] / 60.0
        return df

    def _remove_adduct_mass(self):
        def mass_of_adducts(adduct):
            full_change = 0.0
            split = re.split('(\+|-)', adduct)[1:]
            plus = None
            for anno in split:
                if anno == "+":
                    plus = True
                elif anno == '-':
                    plus = False
                else:
                    cf = self._parse_cf(anno)
                    for k, v in cf.items():
                        if plus:
                            full_change += atomic_mass[k] * v
                        elif plus is None:
                            raise Exception("Invalid Adduct used")
                        else:
                            full_change -= atomic_mass[k] * v
            return full_change

        self._id_df['Adduct_Change'] = self._id_df['Adduct'].apply(mass_of_adducts)
        self._id_df['Precursor m/z'] = self._id_df['Precursor m/z'] - self._id_df["Adduct_Change"]

    @staticmethod
    def _expand_to_adducts(df):
        expanded_df = pd.DataFrame(columns=df.columns)
        for row in df.itertuples():
            try:
                temp_row = df.loc[row[0]]
            except:
                temp_row = None
                print("Oops, what happened")

            formmap = L4D_Analyzer_Converter._parse_cf(row.cf)

            default_mz = temp_row['Precursor m/z']

            if settings.H_adduct:
                hydrogen_adduct = deepcopy(formmap)
                hydrogen_adduct['H'] += 1
                hydrogen_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in hydrogen_adduct.items())
                hydrogen_mz = default_mz + atomic_mass["H"]
                temp_row['Adduct'] = 'M+H'
                temp_row['Adduct_cf'] = hydrogen_cf
                temp_row['Precursor m/z'] = hydrogen_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.Na_adduct:
                sodium_adduct = deepcopy(formmap)
                sodium_adduct['Na'] += 1
                sodium_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in sodium_adduct.items())
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
                    hydrogen_minus_water_cf = ''.join(
                        '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in
                        hydrogen_minus_water.items())
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
                    sodium_minus_water_cf = ''.join(
                        '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in
                        sodium_minus_water.items())
                    sodium_minus_water_mz = default_mz + atomic_mass["Na"] - atomic_mass["H"] * 2 - atomic_mass["O"]
                    temp_row['Adduct'] = 'M+Na-[H2O]'
                    temp_row['Adduct_cf'] = sodium_minus_water_cf
                    temp_row['Precursor m/z'] = sodium_minus_water_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.NH4_adduct:
                ammonium_adduct = deepcopy(formmap)
                ammonium_adduct['N'] += 1
                ammonium_adduct['H'] += 4
                ammonium_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in ammonium_adduct.items())
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
                    ammonium_minus_water_cf = ''.join(
                        '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in
                        ammonium_minus_water.items())
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
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in formate_adduct.items())
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
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in acetate_adduct.items())
                acetate_mz = default_mz + atomic_mass["H"] + atomic_mass["C"] + atomic_mass["O"] * 2 + ELECTRON_MASS
                temp_row['Adduct'] = '+C2H3O2-'
                temp_row['Adduct_cf'] = acetate_cf
                temp_row['Precursor m/z'] = acetate_mz
                expanded_df = expanded_df.append(temp_row, ignore_index=True)

            if settings.e_adduct:
                # electron:
                electron_adduct = deepcopy(formmap)
                electron_cf = ''.join(
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in electron_adduct.items())
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
                    '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in hydrogen_abduct.items())
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
                        '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in water_abduct.items())
                    water_abduct_mz = default_mz - atomic_mass["H"] * 2 - atomic_mass["O"] + ELECTRON_MASS
                    temp_row['Adduct'] = '-H2O'
                    temp_row['Adduct_cf'] = water_abduct_cf
                    temp_row['Precursor m/z'] = water_abduct_mz
                    expanded_df = expanded_df.append(temp_row, ignore_index=True)

        return expanded_df

    @staticmethod
    def _parse_cf(cf):
        import re
        input_string = cf
        formmap = {"C": 0, "H": 0, 'N': 0, "O": 0, "P": 0}
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
    def _write_cf(cf):
        return ''.join(
            '%s%s' % (k, v) if v > 1 else '%s' % k if v == 1 else '' for k, v in cf.items())

    def _set_n_peaks(self):
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
                      sub_list[id_charge_index] * L4D_Analyzer_Converter.PROTON_MASS
            for z in range(settings.min_charge_state, settings.max_charge_state + 1):
                quick_list = copy(sub_list)
                quick_list[precursor_mz_index] = (premass + float(z) * L4D_Analyzer_Converter.PROTON_MASS) / float(z)
                quick_list[id_charge_index] = z
                all_data.append(quick_list)
        self._id_df = pd.DataFrame(all_data, columns=all_columns)

    def _finalize(self):
        header_order = L4D_Analyzer_Converter.correct_header_order

        self._id_df.sort_index(inplace=True)
        self._id_df['cf'] = self._id_df['cf'].apply(lambda x: self._write_cf(self._parse_cf(x)))
        self._id_df['literature_n'] = np.nan
        self._id_df['HMP'] = np.nan
        self._id_df['LMP'] = np.nan
        self._id_df['Lipid Unique Identifier'] = self._id_df['Lipid Name'] + "_" + self._id_df[
            'Precursor Retention Time (sec)'].apply(lambda x: str("{:.3f}".format(x / 60))) + str(self._id_df[
                                                                                                      'Precursor Drift Time (ms)'])
        self._id_df = self._id_df[header_order].rename(
            columns=L4D_Analyzer_Converter.correct_header_names
        )
        self._id_df.sort_values(inplace=True, by=["Precursor Retention Time (sec)", "Precursor Drift Time (ms)"])


def main():
    # file = ["/media/JCPriceLab/BUBBLES/Chad Q/Lipid COnstituent Summation/results_10ppm_w_cf_Aug_9.tsv"]
    file = [tk_get_single_file()]
    converter = L4D_Analyzer_Converter(file)
    converter.convert()
    print("done")

    converter.write(file[0][:-4] + "_ID_File.csv")


if __name__ == "__main__":
    main()
