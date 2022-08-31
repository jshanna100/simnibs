from scipy.io import loadmat
import os
import pandas as pd
import re

version = 4
root_dir = "/home/hannaj/"
data_dir = os.path.join(root_dir, f"simnibs/{version}_results")

res_dirs = next(os.walk(data_dir))[1]

failed_dirs = []
df_vars = ["Project", "ROI", "Condition", "Subject", "Radii", "Mags", "Focs",
           "BestRadius", "BestMag", "BestFoc", "Phi"]
df_dict = {x:[] for x in df_vars}
for rd in res_dirs:
    files = os.listdir(os.path.join(data_dir, rd))
    if "summary_metrics.mat" not in files:
        failed_dirs.append(rd)
        print(f"{rd} does not have final metrics.")
        continue
    re_rd = re.match("(.*)_(.*)__(\d*)__(.*)", rd)
    project, roi, subj, cond = re_rd.groups()
    mat = loadmat(os.path.join(data_dir, rd, "summary_metrics.mat"))

    df_dict["Project"].append(project)
    df_dict["ROI"].append(roi)
    df_dict["Subject"].append(subj)
    df_dict["Condition"].append(cond)
    df_dict["Radii"].append(mat["radius_surround"][0])
    df_dict["Mags"].append(mat["roi_median_r"].flatten())
    df_dict["Focs"].append(mat["focality_r"].flatten())
    df_dict["BestRadius"].append(mat["best_radius"][0][0])
    df_dict["BestMag"].append(mat["roi_median_f"][0][0])
    df_dict["BestFoc"].append(mat["focality_f"][0][0])
    df_dict["Phi"].append(mat["best_phi"][0][0])

df = pd.DataFrame.from_dict(df_dict)
df.to_pickle(os.path.join(data_dir, f"df_{version}.pickle"))
