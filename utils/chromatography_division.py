import pandas as pd
import numpy as np

try:
    from obs.peak import Peak
    from obs.id import ID
    from obs.envelope import Envelope
    from obs.molecule import Molecule
    import deuterater.settings as settings
except: # noqa
    print("Your directories are in a weird place...")
#     from DeuteRater.obs.peak import Peak
#     from DeuteRater.obs.id import ID
#     from DeuteRater.obs.envelope import Envelope
#     from DeuteRater.obs.molecule import Molecule
#     import DeuteRater.deuterater.settings as settings
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
                # BD: Issue with mp.cpu_count() finding too many cores available
                self._n_processors = round(mp.cpu_count() * 0.75)
                # self._n_processors = mp.cpu_count()
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
        parsed[0] = parsed[0].replace("[", "")
        parsed[-1] = parsed[-1].replace("]", "")

        return np.array(parsed)

    @staticmethod
    def divide_molecule(molecule):
        molecule.analyze_peaks()
        molecule.choose_peak()
        molecule.extract_with_chosen_peak()
        return molecule

    @staticmethod
    def handle_molecule(group, settings_path):
        settings.load(settings_path)
        molecule = Molecule(group[0])
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
        for series in df.iterrows():
            row = series[1]
            mzs = ChromatographyDivider.parse_2d_list(row.mzs_list)
            abs = ChromatographyDivider.parse_2d_list(row.intensities_list)
            rts = ChromatographyDivider.parse_1d_list(row.rt_list)
            baselines = ChromatographyDivider.parse_1d_list(row.baseline_list)

            # Remove any data from memory so that I don't kill the memory.
            row.mzs_list = np.nan
            row.intensities_list = np.nan
            row.rt_list = np.nan
            row.baseline_list = np.nan
            df.loc[df.name_check == row.name_check, "mzs_list"] = np.nan
            df.loc[df.name_check == row.name_check, "intensities_list"] = np.nan
            df.loc[df.name_check == row.name_check, "rt_list"] = np.nan
            df.loc[df.name_check == row.name_check, "baseline_list"] = np.nan

            envelopes = list()
            should_plot = False
            num_normal_peaks = abs.shape[1] - settings.peak_lookback - settings.peak_lookahead
            for envelope_num in range(abs.shape[0]):
                try:
                    peaks = [Peak(mzs[envelope_num][i], abs[envelope_num][i], i - 1) for i in range(abs.shape[1])]
                except:
                    peaks = None
                    print("hello")
                envelope = Envelope(peaks, rts[envelope_num], settings.peak_lookback, settings.peak_lookahead)
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
            id._envelopes = envelopes
            if should_plot:
                id.divide_chromatography()

            # Separate a molecule by its Adduct, charge, and mzml_file
            column_list = row.index.isin(["Adduct", "z", "mzml_name"])
            column_list = row.index[column_list]
            column_list = list(column_list)
            column_list.sort()
            unique_identifier = "_".join(row.loc[column_list].astype(str))
            molecule.add_id(id, unique_identifier)

        molecule = ChromatographyDivider.divide_molecule(molecule)
        group_df = molecule.update_output_file(df)

        return group_df

    def divide(self):
        if self.biomolecule_type == "Lipid":
            col_names = ["Lipid Unique Identifier", "rt_list", "mzml_path", "Adduct", "z"]
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
            molecule_groups = df.groupby(by=molecule_group_name)
            molecules = list()

            num_groups = molecule_groups.ngroups
            if settings.debug_level == 0:
                func = partial(self.handle_molecule, settings_path=self.settings_path)
                molecules = self._mp_pool.map(func,
                                              tqdm(molecule_groups, desc="dividing chromatography: ",
                                                   total=num_groups))

                with open("logs.txt", 'a') as log:
                    log.write("Settings Path: " + str(self.settings_path) + "\n")
                    log.write("Input Paths: " + str(self.input_paths) + "\n")
                    log.write("Output Paths: " + str(self.out_paths) + "\n")
                    log.write("Interfile\n")
                    log.write("Biomolecule type: " + self.biomolecule_type + "\n")
                    log.write("DF Columns: " + str(df.columns) + "\n")
                    log.write("Column Names: " + str(col_names) + "\n")
                    log.write("Molecule Group Name: " + str(molecule_group_name) + "\n")
                    log.write("Molecule Groups: " + str(molecule_groups.groups) + "\n")
                    log.write("Molecules: " + str(molecules) + "\n\n")

            elif settings.debug_level >= 1:
                func = partial(self.handle_molecule, settings_path=self.settings_path)
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
    print("This is a utility file. ",
          "If you would like to use the Chromatography Division by itself, ",
          "run the file chromatography_division_sep.py inside the main DeuteRater folder.")


if __name__ == "__main__":
    main()
