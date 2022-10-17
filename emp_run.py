from simnibs import __version__
from simu_func import emp_montage
from itertools import product
import os
from joblib import Parallel, delayed

def build_subject_paths(root_dir, version):
    subjectpaths = []
    subj_dirs = next(os.walk(root_dir))[1]
    for subj_dir in subj_dirs:
        if round(version) == 3:
            m2m_str = f"m2m_{subj_dir}"
        elif round(version) == 4:
            m2m_str = f"m2m_T1w.nii"
        subjectpaths.append(os.path.join(root_dir, subj_dir, m2m_str))
    subj_dict = [{subject_dir:subj_path} for
                 subj_path, subject_dir in zip(subjectpaths, subj_dirs)]
    return subj_dict

version = int(__version__[0])

projs = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
masks = ["P1_rTP-RH", "P2_lPCC-new-LH",  "P3_lTP-LH", "P4_lIFG-LH", "P5_lM1-LH",
         None, "P7_rDLPFCnew-RH", "P8_lDLPFC-LH"]
hemis = ['rh', "lh", "lh", "lh", "lh", None, "rh", "lh"]

# projs = ["P6"]
# masks = [None]
# hemi = [None]

proj_dicts = [{proj:[mask, hemi]} for proj, mask, hemi in zip(projs, masks,
                                                              hemis)]
root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
#root_dir = "/home/jev/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))
subj_dicts = build_subject_paths(data_dir, version)
n_jobs = 6
queue = list(product(subj_dicts, proj_dicts))
#queue = [queue[0]]
Parallel(n_jobs=n_jobs)(delayed(emp_montage)(*q, root_dir) for q in queue)
