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

version = 3
root_dir = "/home/hannaj/"
root_dir = "/home/jev/"
data_dir = join(root_dir, f"simnibs/{version}_results")
fig_dir = join(root_dir, "simnibs/figures")

df = pd.read_pickle(join(data_dir, f"df_{version}.pickle"))

projects = np.sort(df["Project"].unique())

for project in projects:
    this_df = df.query(f"Project=='{project}'")
    subjs = np.sort(this_df["Subject"].unique())

    fig = plt.figure(figsize=(38.4, 21.6))
    plt.suptitle(f"{project} (v{version})", fontsize=28)
    col_n = 10
    if version == 4:
        row_n = 4
        valid_rows = 2
        alpha = 0.5
    else:
        row_n = 5
        valid_rows = 1
        alpha = 1.

    outer = gridspec.GridSpec(row_n, col_n)
    for subj_idx, subj in enumerate(subjs):
        subj_df = this_df.query(f"Subject=='{subj}'")
        if len(subj_df) != valid_rows:
            continue
        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                                                 subplot_spec=outer[subj_idx])

        closest_radii = subj_df.query("Condition=='closest'")["Radii"].values[0]
        closest_mags = subj_df.query("Condition=='closest'")["Mags"].values[0]
        closest_focs = subj_df.query("Condition=='closest'")["Focs"].values[0]
        if version == 4:
            optimal_radii = subj_df.query("Condition=='optimal'")["Radii"].values[0]
            optimal_mags = subj_df.query("Condition=='optimal'")["Mags"].values[0]
            optimal_focs = subj_df.query("Condition=='optimal'")["Focs"].values[0]

            if not np.array_equal(closest_radii, optimal_radii):
                raise ValueError("Radii of closest and optimal not the same.")
            else:
                radii = closest_radii
        else:
            radii = closest_radii

        ax = plt.Subplot(fig, inner[0])
        ax.set_title(subj)
        ax.plot(radii, closest_mags, label="closest", alpha=alpha)
        if version == 4:
            ax.plot(radii, optimal_mags, label="optimal", alpha=alpha)
        ax.set_ylim([0.05, 1.])
        ax.set_xticks(radii)
        ax.set_xticklabels([])
        ax.axhline(0.2, linestyle="--", color="black", alpha=.5)
        if outer[subj_idx].colspan[0] == 0:
            ax.set_ylabel("Mag.", fontsize=18, fontweight="bold")
            ax.set_yticks(np.round(np.arange(0.1, 1., .2), 1))
        else:
            ax.set_yticks([])
        fig.add_subplot(ax)

        ax = plt.Subplot(fig, inner[1])
        ax.plot(radii, closest_focs, label="closest", alpha=.5)
        if version == 4:
            ax.plot(radii, optimal_focs, label="optimal", alpha=.5)
        ax.set_ylim([250, 5000])
        if outer[subj_idx].colspan[0] == 0:
            ax.set_ylabel("Foc.", fontsize=18, fontweight="bold")
            ax.set_yticks(np.arange(500, 7500, 2000))
        else:
            ax.set_yticks([])
        ax.set_xticks(radii)
        if outer[subj_idx].rowspan[0] != (row_n - 1):
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Radii", fontsize=18, fontweight="bold")
        fig.add_subplot(ax)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
    plt.savefig(join(fig_dir, f"focmag_{project}_v{version}.png"))
