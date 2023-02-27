from simnibs import __version__
from simnibs import brainsight, localite
import os
import pickle
from simnibs.utils.file_finder import SubjectFiles
from simnibs.utils.matlab_read import read_mat
from simnibs import mesh_io, sim_struct, __version__
import numpy as np

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
root_dir = "/home/jev/temp/"
data_dir = os.path.join(root_dir, str(round(version)))
subj_dicts = build_subject_paths(data_dir)

projs = ["P1", "P2", "P3", "P4", "P5"]

subj_dicts = build_subject_paths(data_dir)
for subj_dict in subj_dicts:
    subname, subpath = list(subj_dict.keys())[0], list(subj_dict.values())[0]
    subject_files = SubjectFiles(subpath=subpath)
    S = read_mat("/home/jev/temp/4_results/P1_rTP-RH__001DO2022__closest/final/simnibs_simulation_20230215-141824.mat")
    msh_file = subject_files.fnamehead
    m = mesh_io.read_msh(subject_files.fnamehead)
    tmslist = S.add_tmslist()
    names = []
    for elec_idx, elec in enumerate(S.poslists[0].electrode):
        coil = tmslist.add_position()
        coil.centre = elec.centre
        coil.name = elec.name
        if elec.name:
            names.append(elec.name)
        else:
            names.append(f"{elec_idx}")
        coil.pos_ydir = np.array([0,0,0])
    matsimnibs = [p.calc_matsimnibs(msh_file, cap=S.eeg_cap) for p in tmslist.pos]
    brainsight().write(tmslist, f"/home/jev/temp/{subname}_test_lps",
                       names=names, overwrite=True,
                       out_coord_space="World")
    localite().write(tmslist, f"/home/jev/temp/{subname}_test_lps_localite",
                     names=names, overwrite=True,
                     out_coord_space="LPS")

    brainsight().write(tmslist, f"/home/jev/temp/{subname}_test_ras",
                       names=names, overwrite=True,
                       out_coord_space="NifTI:Scanner")
    localite().write(tmslist, f"/home/jev/temp/{subname}_test_ras_localite",
                     names=names, overwrite=True,
                     out_coord_space="RAS")
