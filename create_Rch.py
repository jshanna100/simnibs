import os
import pandas as pd

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

version = 4

root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
root_dir = "/home/jev/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))

subj_dicts = build_subject_paths(data_dir, version)

for subj_dict in subj_dicts:
    try:
        subj_dir = list(subj_dict.values())[0]
        this_path = os.path.join(subj_dir, "eeg_positions",
                                 "EEGcap_incl_cheek_buci_2.csv")
        df = pd.read_csv(this_path, header=None)

        chk_row = df[df[4]=="Lch"].copy()
        flip_x = -1 * chk_row[1].values[0]
        chk_row[1] = flip_x
        chk_row[4] = "Rch"
        # move down a few cm to avoid eye
        chk_row[3] -= 5.
        df = df.append(chk_row)

        chk_row = df[df[4]=="yLch"].copy()
        flip_x = -1 * chk_row[1].values[0]
        chk_row[1] = flip_x
        chk_row[4] = "yRch"
        # move down a few cm to avoid eye
        chk_row[3] -= 5.
        df = df.append(chk_row)

        out_path = os.path.join(subj_dir, "eeg_positions",
                                 "EEGcap_incl_cheek_buci_4.csv")
        df.to_csv(out_path, header=False, index=False)
    except:
        print(f"\n\nCould not add channels for subject {subj_dir}\n\n")
