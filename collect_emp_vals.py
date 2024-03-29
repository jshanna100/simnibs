from scipy.io import loadmat
import os
import pandas as pd
import re

version = 4
root_dir = "/home/hannaj/"
#root_dir = "/home/jev/"
data_dir = os.path.join(root_dir, f"simnibs/{version}_emp")

res_dirs = next(os.walk(data_dir))[1]

failed_dirs = []
df_vars = ["Project", "Subject", "Mag", "Mag", "Foc",
           "Foc", "Summary"]
df_dict = {x:[] for x in df_vars}
for rd in res_dirs:
    files = os.listdir(os.path.join(data_dir, rd))
    if "summary_metrics.mat" not in files:
        failed_dirs.append(rd)
        print(f"{rd} does not have final metrics.")
        continue
    re_rd = re.match("(.*)_(.*)", rd)
    subj, proj = re_rd.groups()
    mat = loadmat(os.path.join(data_dir, rd, "summary_metrics.mat"))

    df_dict["Project"].append(proj)
    df_dict["Subject"].append(subj)
    df_dict["Summary"].append("ROI_Median")
    df_dict["Mag"].append(mat["median"][0][0])
    df_dict["Foc"].append(mat["focality_med"][0][0])
    df_dict["Project"].append(proj)
    df_dict["Subject"].append(subj)
    df_dict["Summary"].append("ROI_Mean")
    df_dict["Mag"].append(mat["mean"][0][0])
    df_dict["Foc"].append(mat["focality_mean"][0][0])

df = pd.DataFrame.from_dict(df_dict)
df.to_pickle(os.path.join(data_dir, f"df_emp_{version}.pickle"))
