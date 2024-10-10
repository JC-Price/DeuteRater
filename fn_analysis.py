import pandas as pd
import tkinter as tk
import tkinter.filedialog as fd
import os


def main():
    root = tk.Tk()
    files = fd.askopenfilenames(parent=root, title='Choose file(s)')

    for file in files:
        df = pd.read_csv(file, sep='\t')
        df = df.loc[df['enrichment'] == 0.05]
        df = df[pd.to_numeric(df['cfn'], errors='coerce').notnull()]
        df['cfn'] = df['cfn'].astype(float)
        avg_std = df['frac_new_combined_std_dev'].mean()
        cfn = df['cfn'].mean()
        basename = os.path.basename(file)

        print(f"{basename} combined average standard deviation: {avg_std}")
        print(f"{basename} combined average fraction new: {cfn}")


if __name__ == "__main__":
    main()
