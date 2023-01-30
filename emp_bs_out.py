from simnibs import __version__
from simnibs import brainsight
from simu_func import emp_montage
from itertools import product
import os
import pickle
from emp_chandefs import prepare_emp

def build_subject_paths(root_dir):
    subjectpaths = []
    subj_dirs = next(os.walk(root_dir))[1]
    for subj_dir in subj_dirs:
        subjectpaths.append(os.path.join(root_dir, subj_dir, f"m2m_{subj_dir}"))
    subj_dict = [{subject_dir:subj_path} for
                 subj_path, subject_dir in zip(subjectpaths, subj_dirs)]
    return subj_dict

version = int(__version__[0])
#root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
root_dir = "/home/jev/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))
subj_dicts = build_subject_paths(data_dir)

projs = ["P1", "P2", "P2-5050", "P3", "P4", "P5", "P5-2ma", "P5-2ma-5050",
         "P6", "P6-3030", "P7", "P8"]

subj_dicts = build_subject_paths(data_dir)
for subj_dict in [subj_dicts[0]]:
    subname, subpath = list(subj_dict.keys())[0], list(subj_dict.values())[0]
    for project in projs:
        S = prepare_emp(project, tms=True)
        S.subpath = subpath
        if "P2" in project or "P6" in project:
            S.eeg_cap = S.subpath + '/eeg_positions' + '/EEGcap_incl_cheek_buci_2.csv'
        breakpoint()
