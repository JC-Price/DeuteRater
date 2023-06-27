import pandas as pd
import numpy as np
import re

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


def correct_mz(df):
    for idx, row in df.iterrows():
        adduct = row["ConsistentAdducts"]
        adduct = re.split(";|\+", adduct)[1]
        if adduct == "H":
            row["mz"] -= atomic_mass["H"]
        elif adduct == "Na":
            row["mz"] -= atomic_mass["Na"]
        elif adduct == "NH4":
            row["mz"] -= (atomic_mass["H"] * 4) + atomic_mass["N"]
        elif adduct == "-H":
            row["mz"] += atomic_mass["H"]
        elif adduct == "CHO2":
            row["mz"] -= atomic_mass["C"] + atomic_mass["H"] + (atomic_mass["O"] * 2)
        elif adduct == "C2H3O2":
            row["mz"] -= (atomic_mass["C"] * 2) + (atomic_mass["H"] * 3) + (atomic_mass["O"] * 2)
        else:
            print(f"CAN'T HANDLE THE ADDUCT {adduct}!!!")
            exit(-1)

        df.at[idx, "mz"] = row["mz"]
    return df


filename = "results-test.csv"
df = pd.read_csv(filename)

results = df.dropna(subset=["ConsistentAdducts"])
results = correct_mz(results)

pre_df = pd.DataFrame()
pre_df["Name"] = results["lipidMolecularSpecies"]
pre_df["RT"] = results["rt"]
pre_df["Mass"] = results["mz"]
pre_df["CCS"] = results['ccs']
try:
    pre_df["DT"] = results['DT']
    pre_df["cf"] = results['Chemical Formula']
except:  # noqa
    print("Need DT and Formula added to file before preprocessing")

pre_df.to_csv("pre_Id_file.csv")
print('hi')

