import numpy as np
import os
import shutil
from scipy.io import savemat
from itertools import product
from joblib import Parallel, delayed
from simnibs import __version__, sim_struct, mesh_io
from simnibs.utils.file_finder import SubjectFiles
import Nx1_stuff
from simu_func import rad_only

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

masks = ["P1_rTP-RH", "P2_lPCC-new-LH",  "P3_lTP-LH", "P4_lIFG-LH", "P5_lM1-LH",
         "P7_rDLPFCnew-RH", "P8_lDLPFC-LH"]
phis = [35., 90., 90., 75., 90., 30., 75.]
hemis = ['rh', "lh", "lh", "lh", "lh", "rh", "lh"]
# mask_dicts contains phis and hemisphere for each ROI
mask_dicts = [{mask:[phi, hemi]} for mask, phi, hemi in zip(masks, phis, hemis)]
conditions = ["closest", "optimal"]

version = int(__version__[0])

root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
#root_dir = "/home/jev/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))
subj_dicts = build_subject_paths(data_dir, version)
n_jobs = 6


print(f"\nVersion {version}\n")

radius_surround = list(np.arange(30, 100, 10))
# properties of centre electrode
EL_center = sim_struct.ELECTRODE()
EL_center.shape = 'ellipse'  # round shape
EL_center.dimensions = [20, 20]  # 30 mm diameter
EL_center.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer

# properties of surround electrodes
EL_surround = sim_struct.ELECTRODE()
EL_surround.shape = 'ellipse'  # round shape
EL_surround.dimensions = [20, 20]  # 30 mm diameter
EL_surround.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer

queue = list(product(subj_dicts, mask_dicts, conditions)) # all combinations of subjects/masks/condition
args = [radius_surround, EL_center, EL_surround, root_dir]
Parallel(n_jobs=n_jobs)(delayed(rad_only)(*q, *args) for q in queue[2:])
