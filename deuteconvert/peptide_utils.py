# TODO: Documentation is great.


def calc_cf(aa_counts, aa_composition_df,z):
    elem_dict = {elem: 0 for elem in aa_composition_df.loc[:, 'C':].columns}
    for aa, count in aa_counts.items():
        for elem in elem_dict.keys():
            elem_dict[elem] += count * aa_composition_df.at[aa, elem]
    elem_dict = positive_cf_change(elem_dict, z)
    return elem_dict

#the calc_cf function calculates elemental composition well
#however it assumes all AAs have pepitde bonds on both ends  
#(which is not true of terminal AAs) and does not include the charge values
#so h2O and one H per charge must be added to get the most accurate cf
#if we go to negative mode or use different adducts, adjust here

#adds 1 O and z+2 H values to the elemental dict.
#takes elemental_dict (dictionary) and z (an int)
#returns the updated elemental dictionary
def positive_cf_change(elemental_dict, z):
    # since we are assuming positve and +H adduct we can just be direct
    elemental_dict["O"] +=1
    elemental_dict["H"] += 2 + z
    return elemental_dict


def calc_theory_mass(cf, elements_df):
    mass = 0
    for elem, count in cf.items():
        mass += count * elements_df.at[elem, 'relative_atomic_mass']
    return mass


def calc_add_n(aa_counts, aa_labeling_dict):
    add_n = 0
    for aa, count in aa_counts.items():
        add_n += count * aa_labeling_dict[aa]
    return add_n
