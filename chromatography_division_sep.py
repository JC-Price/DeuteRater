import os.path

from utils.chromatography_division import ChromatographyDivider


def main():
    import sys
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("error", category=RuntimeWarning)

    def tk_get_files(extension='*', prompt="Select file"):
        from tkinter import filedialog
        from tkinter import Tk
        import os
        root = Tk()
        root.withdraw()
        if extension == '*':
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

    settings_path = os.path.join(os.getcwd(), "resources/peptide_settings.yaml")

    input_files = tk_get_files()
    if len(input_files[0]) == 1:
        input_files = [input_files]
    out_files = [f[:-4] + "_divided.tsv" for f in input_files]

    divider = ChromatographyDivider(settings_path=settings_path,
                                    input_paths=input_files,
                                    out_paths=out_files,
                                    biomolecule_type="Peptide"
                                    )
    divider.divide()


if __name__ == "__main__":
    main()
