import pandas as pd
import numpy as np

no_CD_file = "F:\\DR Testing\\FullData, no CD, short guide file\\frac_new_output.tsv"
CD_file = "F:\\DR Testing\\FullData, interfile, short guide file, mzML change\\frac_new_output.tsv"
intra_file = "F:\\DR Testing\\FullData, intrafile, short guide file, mzML change\\frac_new_output.tsv"

lipid_no_CD_file = "F:\\DR Testing\\LipidData, no CD, short guide file\\frac_new_output.tsv"
lipid_CD_file = "F:\\DR Testing\\LipidData, CD - interfile, short guide file\\frac_new_output.tsv"
lipid_intra_file = "F:\\DR Testing\\LipidData, CD - intrafile, short guide file\\frac_new_output.tsv"

lipid_set_enrich_file = "F:\\DR Testing\\Enrichment test\\frac_new_output.tsv"


def main():
    no_cd = pd.read_csv(no_CD_file, sep="\t")
    cd = pd.read_csv(CD_file, sep="\t")
    intra = pd.read_csv(intra_file, sep='\t')
    
    lipid_no_cd = pd.read_csv(lipid_no_CD_file, sep="\t")
    lipid_cd = pd.read_csv(lipid_CD_file, sep="\t")
    lipid_intra = pd.read_csv(lipid_intra_file, sep='\t')
    lipid_enrich = pd.read_csv(lipid_set_enrich_file, sep='\t')
    
    no_cd = no_cd.loc[no_cd['enrichment'] == 0.0]
    cd = cd.loc[cd['enrichment'] == 0.0]
    intra = intra.loc[intra['enrichment'] == 0.0]
    
    lipid_no_cd = lipid_no_cd.loc[lipid_no_cd['enrichment'] < 0.01]
    lipid_cd = lipid_cd.loc[lipid_cd['enrichment'] < 0.01]
    lipid_intra = lipid_intra.loc[lipid_intra['enrichment'] < 0.01]
    lipid_enrich = lipid_enrich.loc[lipid_enrich['enrichment'] < 0.0501]
    
    # no_cd.to_csv("D:\\DR Testing\\analysis\\frac_new_output_no_chromatography.tsv", sep="\t", index=False)
    # cd.to_csv("D:\\DR Testing\\analysis\\frac_new_output_with_chromatography_inter.tsv", sep="\t", index=False)
    # intra.to_csv("D:\\DR Testing\\analysis\\frac_new_output_with_chromatography_intra.tsv", sep="\t", index=False)

    # lipid_no_cd.to_csv("F:\\DR Testing\\lipid_frac_new_output_no_chromatography.tsv", sep="\t", index=False)
    # lipid_cd.to_csv("F:\\DR Testing\\lipid_frac_new_output_with_chromatography_inter.tsv", sep="\t", index=False)
    # lipid_intra.to_csv("F:\\DR Testing\\lipid_frac_new_output_with_chromatography_intra.tsv", sep="\t", index=False)
    lipid_enrich.to_csv("F:\\DR Testing\\lipid_set_enrichment.tsv", sep="\t", index=False)

    no_cd = no_cd[pd.to_numeric(no_cd['frac_new_combined_std_dev'], errors='coerce').notnull()]
    cd = cd[pd.to_numeric(cd['frac_new_combined_std_dev'], errors='coerce').notnull()]
    intra = intra[pd.to_numeric(intra['frac_new_combined_std_dev'], errors='coerce').notnull()]
    
    lipid_no_cd = lipid_no_cd[pd.to_numeric(lipid_no_cd['frac_new_combined_std_dev'], errors='coerce').notnull()]
    lipid_cd = lipid_cd[pd.to_numeric(lipid_cd['frac_new_combined_std_dev'], errors='coerce').notnull()]
    lipid_intra = lipid_intra[pd.to_numeric(lipid_intra['frac_new_combined_std_dev'], errors='coerce').notnull()]
    lipid_enrich = lipid_enrich[pd.to_numeric(lipid_enrich['frac_new_combined_std_dev'], errors='coerce').notnull()]
    
    no_cd = no_cd[pd.to_numeric(no_cd['cfn'], errors='coerce').notnull()]
    cd = cd[pd.to_numeric(cd['cfn'], errors='coerce').notnull()]
    intra = intra[pd.to_numeric(intra['cfn'], errors='coerce').notnull()]
    
    lipid_no_cd = lipid_no_cd[pd.to_numeric(lipid_no_cd['cfn'], errors='coerce').notnull()]
    lipid_cd = lipid_cd[pd.to_numeric(lipid_cd['cfn'], errors='coerce').notnull()]
    lipid_intra = lipid_intra[pd.to_numeric(lipid_intra['cfn'], errors='coerce').notnull()]
    lipid_enrich = lipid_enrich[pd.to_numeric(lipid_enrich['cfn'], errors='coerce').notnull()]
    
    no_cd['cfn'] = no_cd['cfn'].astype(float)
    cd['cfn'] = cd['cfn'].astype(float)
    intra['cfn'] = intra['cfn'].astype(float)
    
    lipid_no_cd['cfn'] = lipid_no_cd['cfn'].astype(float)
    lipid_cd['cfn'] = lipid_cd['cfn'].astype(float)
    lipid_intra['cfn'] = lipid_intra['cfn'].astype(float)
    lipid_enrich['cfn'] = lipid_enrich['cfn'].astype(float)

    no_cd_average_std = no_cd['frac_new_combined_std_dev'].mean()
    cd_average_std = cd['frac_new_combined_std_dev'].mean()
    intra_average_std = intra['frac_new_combined_std_dev'].mean()
    no_cd_average_fn = no_cd['cfn'].mean()
    cd_average_fn = cd['cfn'].mean()
    intra_average_fn = intra['cfn'].mean()

    lipid_no_cd_average_std = lipid_no_cd['frac_new_combined_std_dev'].mean()
    lipid_cd_average_std = lipid_cd['frac_new_combined_std_dev'].mean()
    lipid_intra_average_std = lipid_intra['frac_new_combined_std_dev'].mean()
    lipid_no_cd_average_fn = lipid_no_cd['cfn'].mean()
    lipid_cd_average_fn = lipid_cd['cfn'].mean()
    lipid_intra_average_fn = lipid_intra['cfn'].mean()

    lipid_enrich_average_std = lipid_enrich['frac_new_combined_std_dev'].mean()
    lipid_enrich_average_fn = lipid_enrich['cfn'].mean()
    
    print(f"No chromatography division average combined fraction new standard deviation: {no_cd_average_std}")
    print(f"No chromatography division average combined fraction new: {no_cd_average_fn}\n")
    
    print(f"Interfile chromatography division average combined fraction new standard deviation: {cd_average_std}")
    print(f"Interfile chromatography division average combined fraction new: {cd_average_fn}\n")
    
    print(f"Intrafile chromatography division average combined fraction new standard deviation: {intra_average_std}")
    print(f"Intrafile chromatography division average combined fraction new: {intra_average_fn}\n")
    
    print(f"Lipid - No chromatography division average combined fraction new standard deviation: {lipid_no_cd_average_std}")
    print(f"Lipid - No chromatography division average combined fraction new: {lipid_no_cd_average_fn}\n")

    print(f"Lipid - Interfile chromatography division average combined fraction new standard deviation: {lipid_cd_average_std}")
    print(f"Lipid - Interfile chromatography division average combined fraction new: {lipid_cd_average_fn}\n")

    print(f"Lipid - Intrafile chromatography division average combined fraction new standard deviation: {lipid_intra_average_std}")
    print(f"Lipid - Intrafile chromatography division average combined fraction new: {lipid_intra_average_fn}\n")

    print(f"Lipid - Set enrichment, intrafile combined fraction new standard deviation: {lipid_enrich_average_std}")
    print(f"Lipid - Set enrichment, intrafile combined average fraction new: {lipid_enrich_average_fn}")


if __name__ == "__main__":
    main()
    