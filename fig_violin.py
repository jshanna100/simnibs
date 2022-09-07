import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.gridspec as gridspec
import seaborn as sns
import os
import numpy as np
from os.path import join
from scipy.stats import mode

import matplotlib
font = {'weight' : 'bold',
        'size'   : 18}
matplotlib.rc('font', **font)

root_dir = "/home/hannaj/"
root_dir = "/home/jev/"

fig_dir = join(root_dir, "simnibs/figures")

data_dir = join(root_dir, "simnibs/3_results")
df_3 = pd.read_pickle(join(data_dir, "df_3.pickle"))
df_3["Version"] = ["3"] * len(df_3)
data_dir = join(root_dir, "simnibs/4_results")
df_4 = pd.read_pickle(join(data_dir, "df_4.pickle"))
df_4["Version"] = ["4"] * len(df_4)
df = pd.concat([df_3, df_4])

projects = np.sort(df["Project"].unique())

# v4: closest vs optimal
best_radii = {}
best_vals = {"Project":[], "Subject":[], "Condition":[], "Mags":[], "Focs":[]}
this_df = df.query("Version=='4'")
# get best radius for groups
for project in projects:
    best_radii[project] = {}
    proj_df = this_df.query(f"Project=='{project}'")
    for cond in ["closest", "optimal"]:
        cond_df = proj_df.query(f"Condition=='{cond}'")
        #best_radius = vers_df["BestRadius"].mode().values[0]
        rads = np.stack(cond_df["Mags"].values)
        error = np.abs(rads - 0.2)
        rad_idx = mode(error.argmin(axis=1))[0][0]
        radii = cond_df["Radii"].values[0]
        best_radius = radii[rad_idx]
        best_radii[project][cond] = best_radius
        subjs = np.sort(cond_df["Subject"].unique())
        for subj in subjs:
            line = cond_df.query(f"Subject=='{subj}'")
            rad_idx = list(line["Radii"].values[0]).index(best_radius)
            best_vals["Project"].append(project)
            best_vals["Subject"].append(subj)
            best_vals["Condition"].append(cond)
            best_vals["Mags"].append(line["Mags"].values[0][rad_idx])
            best_vals["Focs"].append(line["Focs"].values[0][rad_idx])
best_vals = pd.DataFrame.from_dict(best_vals)

fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=best_vals, x="Project", y="Mags", hue="Condition", ax=ax)
ax.set_title("Magnitude: V4 - Closest vs Optimal", fontsize=24, fontweight="bold")
ax.set_xlabel("Project (closest/optimal best radius)", fontsize=24,
              fontweight="bold")
ax.set_ylabel("Magnitude", fontsize=24, fontweight="bold")
ax.set_xticklabels([f"{k} ({v['closest']}/{v['optimal']})"
                    for k, v in best_radii.items()])
plt.savefig(join(fig_dir, "Mag_ClosestOptimal_violin.png"))


fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=best_vals, x="Project", y="Focs", hue="Condition", ax=ax)
ax.set_title("Focality: V4 - Closest vs Optimal", fontsize=24, fontweight="bold")
ax.set_xlabel("Project (closest/optimal best radius)", fontsize=24,
              fontweight="bold")
ax.set_ylabel("Focality", fontsize=24, fontweight="bold")
ax.set_xticklabels([f"{k} ({v['closest']}/{v['optimal']})"
                    for k, v in best_radii.items()])
plt.savefig(join(fig_dir, "Foc_ClosestOptimal_violin.png"))

# closest: 3 vs 4
best_radii = {}
best_vals = {"Project":[], "Subject":[], "Version":[], "Mags":[], "Focs":[]}
this_df = df.query("Condition=='closest'")
# get best radius for groups
for project in projects:
    best_radii[project] = {}
    proj_df = this_df.query(f"Project=='{project}'")
    for version in [3, 4]:
        vers_df = proj_df.query(f"Version=='{version}'")
        #best_radius = vers_df["BestRadius"].mode().values[0]
        rads = np.stack(vers_df["Mags"].values)
        error = np.abs(rads - 0.2)
        rad_idx = mode(error.argmin(axis=1))[0][0]
        radii = vers_df["Radii"].values[0]
        best_radius = radii[rad_idx]
        best_radii[project][version] = best_radius
        subjs = np.sort(vers_df["Subject"].unique())
        for subj in subjs:
            line = vers_df.query(f"Subject=='{subj}'")
            rad_idx = list(line["Radii"].values[0]).index(best_radius)
            best_vals["Project"].append(project)
            best_vals["Subject"].append(subj)
            best_vals["Version"].append(version)
            best_vals["Mags"].append(line["Mags"].values[0][rad_idx])
            best_vals["Focs"].append(line["Focs"].values[0][rad_idx])
best_vals = pd.DataFrame.from_dict(best_vals)

fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=best_vals, x="Project", y="Mags", hue="Version", ax=ax)
ax.set_title("Magnitude: Closest - v3 vs. v4", fontsize=24)
ax.set_xlabel("Project (3/4 best radius)", fontsize=24,
              fontweight="bold")
ax.set_ylabel("Magnitude", fontsize=24, fontweight="bold")
ax.set_xticklabels([f"{k} ({v[3]}/{v[4]})"
                    for k, v in best_radii.items()])
plt.savefig(join(fig_dir, "Mag_3vs4_violin.png"))

fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=best_vals, x="Project", y="Focs", hue="Version", ax=ax)
ax.set_title("Focality: Closest - v3 vs. v4", fontsize=24)
ax.set_xlabel("Project (3/4 best radius)", fontsize=24,
              fontweight="bold")
ax.set_ylabel("Focality", fontsize=24, fontweight="bold")
ax.set_xticklabels([f"{k} ({v[3]}/{v[4]})"
                    for k, v in best_radii.items()])
plt.savefig(join(fig_dir, "Foc_3vs4_violin.png"))

df.to_excel(join(fig_dir, "simnibs_3vs4.xlsx"))
