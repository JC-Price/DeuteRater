import pandas as pd
import numpy as np

try:
    from obs.peak import Peak
    from obs.id import ID
    from obs.envelope import Envelope
    from obs.molecule import Molecule
    import deuterater.settings as settings
except:
    from DeuteRater.obs.peak import Peak
    from DeuteRater.obs.id import ID
    from DeuteRater.obs.envelope import Envelope
    from DeuteRater.obs.molecule import Molecule
    import DeuteRater.deuterater.settings as settings
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
        for series in group[1].iterrows():
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
            group[1].loc[group[1].name_check == row.name_check, "mzs_list"] = np.nan
            group[1].loc[group[1].name_check == row.name_check, "intensities_list"] = np.nan
            group[1].loc[group[1].name_check == row.name_check, "rt_list"] = np.nan
            group[1].loc[group[1].name_check == row.name_check, "baseline_list"] = np.nan

            envelopes = list()
            should_plot = False
            num_normal_peaks = abs.shape[1] - settings.peak_lookback - settings.peak_lookahead
            for envelope_num in range(abs.shape[0]):
                try:
                    peaks = [Peak(mzs[envelope_num][i], abs[envelope_num][i], i - 1) for i in range(abs.shape[1])]
                except:
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
        group_df = molecule.update_output_file(group[1])

        return group_df

    def divide(self):
        if self.how_divide == "Interfile":
            dataframes = list()
            for i in range(len(self.input_paths)):
                infile, outfile = self.input_paths[i], self.out_paths[i]
                df = pd.read_csv(infile, sep="\t")
                df = df.dropna(subset=["mzs_list"])
                df["filename"] = outfile
                dataframes.append(df)
            df = pd.concat(dataframes)
            df["mzml_name"] = df["mzml_path"].apply(os.path.basename)

            possible_rows = list(set(df.columns).intersection({"Adduct", "z", "mzml_name"}))
            possible_rows.sort()
            df['name_check'] = df.apply(lambda x: "_".join(x.loc[possible_rows].astype(str)), axis=1)
            if self.biomolecule_type == "Lipid" or self.biomolecule_type == "Lipids":
                molecule_groups = df.groupby(by="Lipid Unique Identifier")
            else:
                molecule_groups = df.groupby(by=df.columns[0])
            df["Extraction_Updated"] = ""
            df["Extraction_Error"] = ""

            molecules = list()

            if settings.debug_level == 0:
                func = partial(self.handle_molecule, settings_path=self.settings_path)
                molecules = self._mp_pool.map(func,
                                              tqdm(molecule_groups, desc="dividing chromatography: ",
                                                   total=len(molecule_groups)))

            elif settings.debug_level >= 1:
                func = partial(self.handle_molecule, settings_path=self.settings_path)
                for group in tqdm(molecule_groups, desc="dividing chromatography: ", total=len(molecule_groups)):
                    molecules.append(func(group))

            df = pd.concat(molecules)

            file_groups = df.groupby(by="filename")
            for file_group in file_groups:
                file_group[1].to_csv(file_group[0], sep="\t", index=False)

        else:
            for i in tqdm(range(len(self.out_paths)), total=len(self.out_paths), desc="dividing files: "):
                infile, outfile = self.input_paths[i], self.out_paths[i]
                df = pd.read_csv(infile, sep='\t')
                df['filename'] = outfile
                df = df.dropna(subset=["mzs_list"])
                possible_rows = list(set(df.columns).intersection({"Adduct", "z", "mzml_name"}))
                possible_rows.sort()
                df['name_check'] = df.apply(lambda x: "_".join(x.loc[possible_rows].astype(str)), axis=1)
                if self.biomolecule_type == "Lipid":
                    molecule_groups = df.groupby(by="Lipid Unique Identifier")
                else:
                    molecule_groups = df.groupby(by=df.columns[0])
                df["Extraction_Updated"] = ""
                df["Extraction_Error"] = ""
                molecules = list()

                if settings.debug_level == 0:
                    func = partial(self.handle_molecule, settings_path=self.settings_path)
                    molecules = self._mp_pool.map(func,
                                                  tqdm(molecule_groups,
                                                       desc="dividing chromatography: ",
                                                       total=len(molecule_groups),
                                                       leave=False))

                elif settings.debug_level >= 1:
                    for group in tqdm(molecule_groups, desc="dividing chromatography: ", total=len(molecule_groups),
                                      leave=False):
                        molecules.append(self.handle_molecule(group))

                df = pd.concat(molecules)

                df.to_csv(outfile, sep='\t', index=False)


def main():
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
        if len(filename) == 1:
            return filename[0]
        return filename  # A string representing the file path

    settings_path = "../resources/settings.yaml"

    input_files = tk_get_files()
    if len(input_files[0]) == 1:
        input_files = [input_files]
    out_files = [f[:-4] + "_divided.tsv" for f in input_files]

    divider = ChromatographyDivider(settings_path=settings_path,
                                    input_paths=input_files,
                                    out_paths=out_files,
                                    biomolecule_type="Lipid"
                                    )
    divider.divide()


if __name__ == "__main__":
    main()
