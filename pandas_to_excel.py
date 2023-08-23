import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.gridspec as gridspec
import seaborn as sns
import os
import numpy as np
from os.path import join

import matplotlib
font = {'weight' : 'bold',
        'size'   : 18}
matplotlib.rc('font', **font)

root_dir = "/home/hannaj/"
root_dir = "/home/jev/"

fig_dir = join(root_dir, "simnibs/figures")

# magfoc
data_dir = join(root_dir, "simnibs/3_results")
df_3 = pd.read_pickle(join(data_dir, "df_3.pickle"))
df_3["Version"] = pd.Series(np.ones(len(df_3))*3)
data_dir = join(root_dir, "simnibs/4_results")
df_4 = pd.read_pickle(join(data_dir, "df_4.pickle"))
df_4["Version"] = pd.Series(np.ones(len(df_4))*4)
df = pd.concat([df_3, df_4])
df = df.query("Condition!='closest_bone'")

mags = np.vstack(df["Mags"].values)
focs = np.vstack(df["Focs"].values)

for idx in range(mags.shape[1]):
    df[f"Mag{idx+1}"] = mags[:, idx]
for idx in range(focs.shape[1]):
    df[f"Foc{idx+1}"] = focs[:, idx]

df.to_excel(join(root_dir, "temp", "radial_magfoc.xlsx"))

# memoslap version
data_dir = join(root_dir, "simnibs")
df = pd.read_pickle(join(data_dir, "df_test.pickle"))

mags_med = np.vstack(df["Mags_med"].values)
mags_sq = np.vstack(df["Mags_sq"].values)
focs = np.vstack(df["Focs"].values)

for idx in range(mags_med.shape[1]):
    df[f"Mags_med{idx+1}"] = mags_med[:, idx]
for idx in range(mags_sq.shape[1]):
    df[f"Mags_sq{idx+1}"] = mags_sq[:, idx]
for idx in range(focs.shape[1]):
    df[f"Foc{idx+1}"] = focs[:, idx]

df.to_excel(join(root_dir, "temp", "radial_magfoc_memo.xlsx"))

# old montages
data_dir = join(root_dir, "simnibs/3_emp")
df_3 = pd.read_pickle(join(data_dir, "df_emp_3.pickle"))
df_3["Version"] = ["3"] * len(df_3)
data_dir = join(root_dir, "simnibs/4_emp")
df_4 = pd.read_pickle(join(data_dir, "df_emp_4.pickle"))
df_4["Version"] = ["4"] * len(df_4)
df = pd.concat([df_3, df_4], ignore_index=True)
df = df.query("Summary=='ROI_Median'")
df.to_excel(join(root_dir, "temp", "old_montage.xlsx"))
