import pandas as pd
import numpy as np  # noqa
from copy import copy
import multiprocessing as mp

from functools import partial

from tqdm import tqdm  # noqa: 401

from utils.emass import parse_cf, emass, normalize

import deuterater.settings as settings

#$at this point we have far too many columns to be practical
#$this is not a problem for theory since that is largely a combination of
#$the table and the extracted files, so it is a nice place to store that data
#$at this point we just need the bare essentials for calculation to ease 
#$troubleshooting and visibility
columns_to_drop = ["Precursor Retention Time (sec)", "rt_start", "rt_end", 
                   "rt_width", "Precursor m/z", "Identification Charge",
                   "id_index", "lookback_mzs", "lookback_abundances",
                   "lookahead_mzs", "lookahead_abundances", "rt_min",
                   "rt_max", "baseline_signal", "mads", 
    ]

peptide_extra_columns_to_drop = ["quality", "avg_ppm", "start_loc", 
                                 "end_loc", "species", 'gene_name',
                                 "protein_existence", "sequence_version"]


#itertuples doesn't like some characters
#these characters include () and space
#this is because the elements are called using periods
#therefore we will declare any particularly problematic.  
#we'll force that to change here 
#$if we're swapping between literature_n and theory_n use this dictionary
protein_itertuple_renamer = {
    "Protein ID": "Protein_ID",
    "Protein Name": "Protein_Name",
    "Homologous Proteins": "Homologous_Proteins",
    "n_isos": "num_peaks",
    "literature_n": "n_value"
    }

lipid_itertuple_renamer = {
    "n_isos": "num_peaks",
    "literature_n": "n_value",
    "Lipid Unique Identifier": "Lipid_Unique_Identifier",
    "Lipid Name": "Lipid_Name"
    }

class FractionNewCalculator():
    def __init__(self, model_path, out_path, settings_path, biomolecule_type):
        settings.load(settings_path)
        self.settings_path = settings_path
        if biomolecule_type == "Peptide":
            itertuple_renamer = copy(protein_itertuple_renamer)
        elif biomolecule_type == "Lipid":
            itertuple_renamer = copy(lipid_itertuple_renamer)
               
        self.biomolecule_type = biomolecule_type
        
        self.error = ""
        
        #$get the number of cores we're using for multiprocessing
        if settings.recognize_available_cores is True:
            self._n_partitions = mp.cpu_count()
        else:
            self._n_partitions = settings.n_partitions
        
        self._mp_pool = mp.Pool(self._n_partitions)
        
        #$ adjust the itertuple renamer to use the proper n_value
        if settings.use_empir_n_value:
            itertuple_renamer['literature_n'] = 'literature_n'
            itertuple_renamer['empir_n'] = 'n_value' 
        self.model = pd.read_csv(model_path, sep='\t')
        self.model.rename(columns = itertuple_renamer, inplace = True)
        self.out_path = out_path

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )

    def generate(self):
        #these need to be removed
        #cf_col = 'cf'
        #nval_col = 'n_value'
        #enrichment_col = 'enrichment'
        
        #$just drop the unneeded columns
        self.model = self.model.drop(columns_to_drop, axis =1, errors = "ignore")
        if self.biomolecule_type == "Peptide":
            self.model = self.model.drop(peptide_extra_columns_to_drop, axis =1,
                                         errors = "ignore")
        
        if self.model["n_value"].isna().all():
            self.error = ("N value is not present.  "
                          "Please ensure the n value has been provided in the "
                          "literature_n column or has been calculated in the "
                          "n_value column. If the proper column is filled "
                          "please check your settings for using the proper "
                          "column")
            return
        
        if settings.use_abundance:
            self.model = self.model.assign(
                theory_unlabeled_abunds="",
                theory_labeled_abunds="",
                normalized_empirical_abundances = "",
                low_labeling_peaks = "",
                frac_new_abunds="",
                frac_new_abunds_std_dev  = "",
                afn = ""
                )
        if settings.use_neutromer_spacing:
            self.model = self.model.assign(
                observed_neutral_masses = "",
                theory_unlabeled_mzs="",
                theory_labeled_mzs="",
                frac_new_mzs="",
                frac_new_mzs_outlier_checked = "",
                frac_new_mzs_std_dev  = "",
                nsfn = ""
                )
        if settings.use_abundance and settings.use_neutromer_spacing:
            self.model = self.model.assign(
                frac_new_combined = "",
                frac_new_combined_outlier_checked = "",
                frac_new_combined_std_dev = "",
                cfn = ""  
                )
        #$now that we have dropped useless columns we can start the multiprocessing
        #$don't use np.split or it will error if chunks are not equal sized
        model_pieces = np.array_split(self.model, self._n_partitions)
        func= partial(FractionNewCalculator._mp_prepare)
        func = partial(func, settings_path = self.settings_path)
        func = partial(func, biomolecule_type = self.biomolecule_type)
        
        results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, model_pieces),
                    total=len(model_pieces)
                )
            
            )
        
        self.model = pd.concat(results)
        
        self._mp_pool.close()
        self._mp_pool.join()
        
    #$generic error output for filters based on n values or other criteria
    #$preferably will do before the calculations so can leave all blank except 
    #$ actual measurement
    @staticmethod
    def _error_method(df, row, error_message):
        if settings.use_abundance:
            df.at[row.Index, "afn"] = error_message
        if settings.use_neutromer_spacing:
            df.at[row.Index, "nsfn"] = error_message
        if settings.use_abundance and settings.use_neutromer_spacing:
            df.at[row.Index, "cfn"] = error_message
        return df

    @staticmethod
    def _mp_prepare(df, settings_path, biomolecule_type):
        settings.load(settings_path)
        #$can start with itertuples.  if need be can swap to apply
        for row in df.itertuples(index=True):
            #$Do any initial filtering to save time calculating
            if row.n_value < settings.min_allowed_n_values:
                df = FractionNewCalculator._error_method(df, row, 
                    "N value is less than {}".format(
                    settings.min_allowed_n_values))
                continue
            if biomolecule_type == "Peptide" and \
                    len(row.Sequence) < settings.min_aa_sequence_length:
                        
                
                df = FractionNewCalculator._error_method(df, row, 
                    "Fewer than {} amino acids".format(
                    settings.min_aa_sequence_length))
                continue
            
            #$if the user chooses 0 as enrichment it will cause divide by zero
            #$error later, so we will use the settings to force it
            #$result should be close to 0 change anyway
            if row.enrichment != 0.0:
                use_enrich = row.enrichment
            else:       
                use_enrich = settings.enrichement_of_zero
                
            num_h, parsed_cf = parse_cf(row.cf)
            _, enriched_results = emass(
                parsed_cf=parsed_cf,
                n_list=[0, row.n_value],
                n_H=num_h,
                low_pct=0,
                high_pct=use_enrich,
                num_peaks=row.num_peaks,
                testing=False
            )
            e_mzs, e_abunds = enriched_results
            
            
            if settings.use_abundance:
                FractionNewCalculator._prepare_row(df, row, "abunds", e_abunds)
                normalized_empirical_abunds = \
                    FractionNewCalculator._normalize_abundances(
                        row.abundances[1:-1].split(", "))
                df.at[row.Index, 'normalized_empirical_abundances']= \
                    normalized_empirical_abunds
            
                theory_abund_deltas = FractionNewCalculator._calculate_deltas(
                    e_abunds.loc[1][1:], e_abunds.loc[0][1:]
                )
                empirical_abund_deltas = FractionNewCalculator._calculate_deltas(
                    [float(x) for x in normalized_empirical_abunds.split(", ")],
                    e_abunds.loc[0][1:]
                )
                theory_abund_deltas, empirical_abund_deltas, removed_peaks = \
                    FractionNewCalculator._trim_abunds(theory_abund_deltas, empirical_abund_deltas)
                
                df.at[row.Index, "low_labeling_peaks"] = removed_peaks
                
                #$don't need to break if we only have one or zero. combined and
                #$spacing are still fine, we just shouldn't do the other calculations
                #$here
                if len(theory_abund_deltas) <2:
                    df.at[row.Index, "afn"] = \
                        f"Insufficient peaks with theory above {settings.minimum_abund_change}"
                    all_frac_new_abunds = []
                else:
                    all_frac_new_abunds = FractionNewCalculator._calculate_fractions(
                        empirical_abund_deltas, theory_abund_deltas
                    )
                    df = FractionNewCalculator.final_calculations(df, row, "abunds", all_frac_new_abunds,
                                             "afn", theory_abund_deltas[0])
            
            if settings.use_neutromer_spacing:     
                FractionNewCalculator._prepare_row(df, row, "mzs", e_mzs)
                theory_mz_deltas = FractionNewCalculator._calculate_deltas(
                   FractionNewCalculator._mz_deltas(e_mzs.loc[1][1:]), 
                   FractionNewCalculator._mz_deltas(e_mzs.loc[0][1:])
                )
                #need to trim the parenthesis from the row.mzs string
                #also turn it into floats to do math on it
                observed_mzs = [float(x) for x in row.mzs[1:-1].split(", ")]
                observed_neutral_mass = [y * int(row.z)
                                         for y in observed_mzs]
                df.at[row.Index, "observed_neutral_masses"] = \
                    ", ".join([str(x) for x in observed_neutral_mass])
                empirical_mz_deltas = FractionNewCalculator._calculate_deltas(
                    FractionNewCalculator._mz_deltas(observed_neutral_mass), 
                    FractionNewCalculator._mz_deltas(e_mzs.loc[0][1:])
                )
                all_frac_new_mzs = FractionNewCalculator._calculate_fractions(
                    empirical_mz_deltas, theory_mz_deltas
                )
                df = FractionNewCalculator.final_calculations(df, row, "mzs", all_frac_new_mzs,
                                         "nsfn")
            
            if settings.use_neutromer_spacing and settings.use_abundance:
               all_frac_new_combined = all_frac_new_abunds + all_frac_new_mzs
               df = FractionNewCalculator.final_calculations(df, row, "combined", all_frac_new_combined,
                                         "cfn")
        return df
            
    @staticmethod
    def _series_to_string(pd_series):
        string_values = [str(a) for a in pd_series.values]
        return " ".join(string_values)
    
    @staticmethod
    def _mz_deltas(A):
        return [mz_value - A[0] for mz_value in A[1:]]

    @staticmethod
    def _calculate_deltas(A, B):
        return [a - b for a, b in zip(A, B)]

    @staticmethod
    def _calculate_fractions(empirical_deltas, theory_deltas):
        # 
        return [a / b for a, b in zip(empirical_deltas, theory_deltas)]
    
    @staticmethod
    def _normalize_abundances(abundance_list):
           normalized_list = normalize([float(a) for a in abundance_list])
           return ", ".join([str(n) for n in normalized_list])
    
   #$need to remove the abundances with a low maximum. usually caused by the
   #$peak rising and falling. since we are calculating based on unidirectional
   #$change this can cause problems for the std dev filter
    @staticmethod 
    def _trim_abunds (theory, empirical):
        drop_list, new_theory, new_empirical = [], [], []
        for t in range(len(theory)):
            if abs(theory[t]) < settings.minimum_abund_change:
                drop_list.append(f'M{t}')
            else:
                new_theory.append(theory[t])
                new_empirical.append(empirical[t])
        return new_theory, new_empirical, ", ".join(drop_list)
    
    #$for right now this is not necessary but we may need to expand it later
    #$ can add a length check or outlier checkto this later
    @staticmethod
    def _std_dev_calculator(list_to_std_dev):
        return np.std(list_to_std_dev, ddof =1)
    
    #$need to add an outlier check to the neutromer spacing (mzs)
    #$as well as any combined data.  This is a median absolute 
    @staticmethod
    def _mad_outlier_check(data_to_check, z_score_cutoff):
        initial_median = np.median(data_to_check) #only calculate once
        differences = [abs(x - initial_median) for x in data_to_check]
        mad = np.median(differences)
        zscores = [.6745 *d/mad for d in differences]
        #$even if we could do a list comprehension for the next step
        #$it would be too long for easy readibility
        good_turnovers = []
        for z in range(len(zscores)):
            if zscores[z] < z_score_cutoff:
                good_turnovers.append(data_to_check[z])
        if len(good_turnovers) > 1:
            return good_turnovers
        else: return data_to_check
        
    #$initial columns used for abundance and spacing.  currently equations are
    #$too different to put in but may do later
    @staticmethod
    def _prepare_row(full_df, row, keyword, df):
        # need to do [1:] in order to cut out the n_D 
        full_df.at[row.Index, 'theory_unlabeled_{}'.format(keyword)] = \
            FractionNewCalculator._series_to_string(df.loc[0][1:])
        full_df.at[row.Index, 'theory_labeled_{}'.format(keyword)] = \
            FractionNewCalculator._series_to_string(df.loc[1][1:])
        return full_df
    
    #$compresses the code to made the last few output columns
    def final_calculations(df, row, keyword, data, final_column,
                           m0_max_delta = 1):
        df.at[row.Index, 'frac_new_{}'.format(keyword)] = \
                ", ".join([str(r) for r in data])
        
        if keyword != "abunds":
            #$ need to do further analysis on outlier check for these
            data = FractionNewCalculator._mad_outlier_check(data, settings.zscore_cutoff)
            df.at[row.Index, 'frac_new_{}_outlier_checked'.format(
                keyword)] = ", ".join([str(r) for r in data])
            df.at[row.Index, 'frac_new_{}_std_dev'.format(keyword)] = \
                str(FractionNewCalculator._std_dev_calculator(data))
            df.at[row.Index, final_column] = str(np.median(data))
        
        else:
            df.at[row.Index, 'frac_new_{}_std_dev'.format(keyword)] = \
                str(FractionNewCalculator._std_dev_calculator(data))
            if abs(m0_max_delta) < settings.min_allowed_abund_max_delta:
                 df.at[row.Index, final_column] = \
                     "Low Theoretical M0 delta possible. Point rejected"
            else:
                df.at[row.Index, final_column] = str(data[0])
        return df