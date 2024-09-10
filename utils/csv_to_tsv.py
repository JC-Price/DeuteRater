import pandas as pd

original_file_path = "F:\\DeuteRater Testing\\Testing Data for Brad\\Proteins\\Full Test Data\\BCM1_NZ1100_ID_file_shortest.csv"
new_file_path = "F:\\DeuteRater Testing\\Testing Data for Brad\\Proteins\\Full Test Data\\BCM1_NZ1100_Test_ID_File.csv"


def main():
    df = pd.read_csv(original_file_path)
    df.to_csv(new_file_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
