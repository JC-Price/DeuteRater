import pandas as pd
import tkinter as tk
import tkinter.filedialog as fd
import os


def main():
    # gather files
    root = tk.Tk()
    files = fd.askopenfilenames(parent=root, title='Choose file(s)')

    for file in files:
        df = pd.read_csv(file, sep='\t')

        # determine if day zero rows had a set enrichment
        day_zero = df.loc[df['time'] < 0.1]
        if day_zero['enrichment'].iloc[0] == 0.05:
            df = df.loc[df['enrichment'] == 0.05]
        else:
            df = df.loc[df['enrichment'] < 0.01]

        # convert cfn to float and ignore cells with text
        df = df[pd.to_numeric(df['cfn'], errors='coerce').notnull()]
        df['cfn'] = df['cfn'].astype(float)

        # calculate average std_dev and fraction new
        avg_std = df['frac_new_combined_std_dev'].mean()
        cfn = df['cfn'].mean()

        # print results
        basename = os.path.basename(file)
        print(f"{basename} combined average standard deviation: {round(avg_std, 4)}")
        print(f"{basename} combined average fraction new: {round(cfn, 4)}\n")


if __name__ == "__main__":
    main()
