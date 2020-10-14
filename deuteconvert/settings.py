import yaml
import traceback
import os

from pathlib import Path


mass_cutoffs: object
rt_proximity_tolerance: int
mz_proximity_tolerance: int
#required_unique: int
start_time: float
study_type: str
aa_elem_comp_path: str
aa_label_path: str
elems_path: str
ptms_path: str
min_charge_state: int
max_charge_state: int

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
resource_location = os.path.join(location, "resources")

# TODO: add quick explanation of how this works, inc. 'global' doc link
def load(settings_path):
    # NOTE: Look at the python documentation for the 'global' statement if you
    #       Want to understand how this module works
    try:
        settings_path = Path(settings_path)
        with settings_path.open('r') as f:
            s = yaml.load(f, Loader=yaml.FullLoader)

        global mass_cutoffs
        mass_cutoffs = sorted(s['mass_cutoffs'])

        global rt_proximity_tolerance
        rt_proximity_tolerance = s['rt_proximity_tolerance']

        global mz_proximity_tolerance
        mz_proximity_tolerance = s['mz_proximity_tolerance']

        global start_time
        start_time = s['start_time']

        global study_type
        study_type = s['study_type']

        global aa_elem_comp_path
        aa_elem_comp_path = os.path.join(resource_location, 
                                         s['aa_elemental_composition_path'])

        global aa_label_path
        aa_label_path = os.path.join(resource_location, 
                                     s['aa_labeling_sites_path'])

        global elems_path
        elems_path = os.path.join(resource_location, s['elements_path'])

        global ptms_path
        ptms_path = os.path.join(resource_location,
                                 s['post_translational_modifications_path'])

        global min_charge_state
        min_charge_state = s['min_charge_state']

        global max_charge_state
        max_charge_state = s['max_charge_state']

    except Exception as e:
        print(e)
        traceback.print_tb(e.__traceback__)


def freeze(path=None):
    settings_dict = {
        'mass_cutoffs': mass_cutoffs,
        'rt_proximity_tolerance': rt_proximity_tolerance,
        'mz_proximity_tolerance': mz_proximity_tolerance,
        'start_time': start_time,
        'study_type': study_type,
        'aa_elemental_composition_path': aa_elem_comp_path,
        'aa_labeling_sites_path': aa_label_path,
        'elements_path': elems_path,
        'post_translational_modifications_path': ptms_path
    }
    if path:
        with open(path, 'w') as frozen_settings_file:
            yaml.dump(
                data=settings_dict,
                stream=frozen_settings_file,
                canonical=False
            )
    else:
        print(yaml.dump(data=settings_dict, canonical=False))
