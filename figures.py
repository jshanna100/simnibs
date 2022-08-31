import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.gridspec as gridspec
import seaborn as sns
import os
import numpy as np

version = 4
root_dir = "/home/hannaj/"
data_dir = os.path.join(root_dir, f"simnibs/{version}_results")

df = pd.read_pickle(os.path.join(data_dir, f"df_{version}.pickle"))

projects = np.sort(df["Project"].unique())

for project in projects:
    this_df = df.query(f"Project=='{project}'")
    subjs = np.sort(this_df["Subject"].unique())

    fig = plt.figure(figsize=(28.4, 19.2))
    plt.suptitle(f"{project} (v{version})", fontsize=28)
    row_n, col_n = 5, 10
    outer = gridspec.GridSpec(row_n, col_n)
    for subj_idx, subj in enumerate(subjs):
        subj_df = this_df.query(f"Subject=='{subj}'")
        if len(subj_df) != 2:
            continue
        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                                                 subplot_spec=outer[subj_idx])

        closest_radii = subj_df.query("Condition=='closest'")["Radii"].values[0]
        closest_mags = subj_df.query("Condition=='closest'")["Mags"].values[0]
        closest_focs = subj_df.query("Condition=='closest'")["Focs"].values[0]
        optimal_radii = subj_df.query("Condition=='optimal'")["Radii"].values[0]
        optimal_mags = subj_df.query("Condition=='optimal'")["Mags"].values[0]
        optimal_focs = subj_df.query("Condition=='optimal'")["Focs"].values[0]

        if not np.array_equal(closest_radii, optimal_radii):
            raise ValueError("Radii of closest and optimal not the same.")
        else:
            radii = closest_radii

        ax = plt.Subplot(fig, inner[0])
        ax.set_title(subj)
        ax.plot(radii, closest_mags, label="closest", alpha=.8)
        ax.plot(radii, optimal_mags, label="optimal", alpha=.8)
        ax.set_ylim([0.05, 0.9])
        ax.set_xticks(radii)
        ax.set_xticklabels([])
        ax.axhline(0.2, linestyle="--", color="black", alpha=.5)
        if outer[subj_idx].colspan[0] == 0:
            ax.set_ylabel("Mag.")
            ax.set_yticks(np.round(np.arange(0.1, .9, .2), 1))
        else:
            ax.set_yticks([])
        fig.add_subplot(ax)

        ax = plt.Subplot(fig, inner[1])
        ax.plot(radii, closest_focs, label="closest", alpha=.8)
        ax.plot(radii, optimal_focs, label="optimal", alpha=.8)
        ax.set_ylim([250, 5000])
        if outer[subj_idx].colspan[0] == 0:
            ax.set_ylabel("Foc")
            ax.set_yticks(np.arange(500, 5000, 1000))
        else:
            ax.set_yticks([])
        ax.set_xticks(radii)
        if outer[subj_idx].rowspan[0] != (row_n - 1):
            ax.set_xticklabels([])
        fig.add_subplot(ax)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
