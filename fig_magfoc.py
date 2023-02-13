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

# best_rads_dict = {"Version":[], "Project":[], "Radius":[]}
#
#
# # compare 3 and 4 with across subject averages
# data_dir = join(root_dir, "simnibs/3_results")
# df_3 = pd.read_pickle(join(data_dir, "df_3.pickle"))
# df_3["Version"] = pd.Series(np.ones(len(df_3))*3)
# data_dir = join(root_dir, "simnibs/4_results")
# df_4 = pd.read_pickle(join(data_dir, "df_4.pickle"))
# df_4["Version"] = pd.Series(np.ones(len(df_4))*4)
# dff = pd.concat([df_3, df_4])
# df = dff.query("Condition=='closest'")
# projects = np.sort(df["Project"].unique())
# mag_fig, mag_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
# plt.suptitle("Magnitude: Closest - 3 vs 4")
# foc_fig, foc_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
# plt.suptitle("Focality: Closest - 3 vs 4")
# mag_axes = [ax for axe in mag_axes for ax in axe]
# foc_axes = [ax for axe in foc_axes for ax in axe]
# for proj_idx, project in enumerate(projects):
#     this_df = df.query(f"Project=='{project}'")
#     # build a new df in a format where this can be easily plotted
#     df_dict = {"Subject":[], "Radius":[], "Magnitude":[], "Focality":[],
#                "Version":[]}
#     radii = this_df.iloc[0]["Radii"]
#     for row_idx, row in this_df.iterrows():
#         for rad_idx, radius in enumerate(radii):
#             df_dict["Subject"].append(row["Subject"])
#             df_dict["Radius"].append(radius)
#             df_dict["Magnitude"].append(row["Mags"][rad_idx])
#             df_dict["Focality"].append(row["Focs"][rad_idx])
#             df_dict["Version"].append(row["Version"])
#     temp_df = pd.DataFrame.from_dict(df_dict)
#
#     sns.lineplot(data=temp_df, x="Radius", y="Magnitude", hue="Version",
#                  ax=mag_axes[proj_idx], palette=["r", "g"])
#     sns.lineplot(data=temp_df, x="Radius", y="Focality", hue="Version",
#                  ax=foc_axes[proj_idx] , palette=["r", "g"])
#     mag_axes[proj_idx].set_title(project)
#     mag_axes[proj_idx].set_ylim([.05, .7])
#     foc_axes[proj_idx].set_ylim([1000, 9000])
#     foc_axes[proj_idx].set_title(project)
#
#     # best radius, 3
#     rad_means = temp_df.query("Version==3.").groupby("Radius")["Magnitude"].mean()
#     best_rad = rad_means[rad_means>0.2].idxmin()
#     mag_axes[proj_idx].plot(best_rad, rad_means[best_rad], marker='o',
#                             color="r", markersize=10)
#     mag_axes[proj_idx].axhline(0.2, color="r", linestyle="--")
#     best_rads_dict["Version"].append("3_closest")
#     best_rads_dict["Project"].append(project)
#     best_rads_dict["Radius"].append(best_rad)
#
#     # best radius, 4 closeset
#     rad_means = temp_df.query("Version==4.").groupby("Radius")["Magnitude"].mean()
#     best_rad = rad_means[rad_means>0.3].idxmin()
#     mag_axes[proj_idx].plot(best_rad, rad_means[best_rad], marker='o',
#                             color="g", markersize=10)
#     mag_axes[proj_idx].axhline(0.3, color="g", linestyle="--")
#     best_rads_dict["Version"].append("4_closest")
#     best_rads_dict["Project"].append(project)
#     best_rads_dict["Radius"].append(best_rad)
#
# # because P6 is missing
# mag_axes[-1].axis("off")
# foc_axes[-1].axis("off")
# mag_fig.savefig(join(fig_dir, "3vs4_mag_projavg.png"))
# foc_fig.savefig(join(fig_dir, "3vs4_foc_projavg.png"))
#
# # compare version 4 closest and optimal with across subject averages
# data_dir = join(root_dir, "simnibs/4_results")
# df = pd.read_pickle(join(data_dir, "df_4.pickle"))
# projects = np.sort(df["Project"].unique())
# mag_fig, mag_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
# plt.suptitle("Magnitude: Version 4 - closest vs optimal")
# foc_fig, foc_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
# plt.suptitle("Focality: Version 4 - closest vs optimal")
# mag_axes = [ax for axe in mag_axes for ax in axe]
# foc_axes = [ax for axe in foc_axes for ax in axe]
# for proj_idx, project in enumerate(projects):
#     this_df = df.query(f"Project=='{project}'")
#     # build a new df in a format where this can be easily plotted
#     df_dict = {"Subject":[], "Radius":[], "Magnitude":[], "Focality":[],
#                "Condition":[]}
#     radii = this_df.iloc[0]["Radii"]
#     for row_idx, row in this_df.iterrows():
#         for rad_idx, radius in enumerate(radii):
#             df_dict["Subject"].append(row["Subject"])
#             df_dict["Radius"].append(radius)
#             df_dict["Magnitude"].append(row["Mags"][rad_idx])
#             df_dict["Focality"].append(row["Focs"][rad_idx])
#             df_dict["Condition"].append(row["Condition"])
#     temp_df = pd.DataFrame.from_dict(df_dict)
#
#     sns.lineplot(data=temp_df, x="Radius", y="Magnitude", hue="Condition",
#                  ax=mag_axes[proj_idx], palette=["r", "g"],
#                  hue_order=["optimal", "closest"])
#     sns.lineplot(data=temp_df, x="Radius", y="Focality", hue="Condition",
#                  ax=foc_axes[proj_idx] , palette=["r", "g"],
#                  hue_order=["optimal", "closest"])
#     mag_axes[proj_idx].set_title(project)
#     mag_axes[proj_idx].set_ylim([.05, .7])
#     foc_axes[proj_idx].set_ylim([1000, 9000])
#     foc_axes[proj_idx].set_title(project)
#
#     # best radius, 4 closest
#     rad_means = temp_df.query("Condition=='closest'").groupby("Radius")["Magnitude"].mean()
#     best_rad = rad_means[rad_means>0.3].idxmin()
#     mag_axes[proj_idx].plot(best_rad, rad_means[best_rad], marker='o',
#                             color="g", markersize=10)
#     mag_axes[proj_idx].axhline(0.3, color="black", linestyle="--")
#     ## already did this above
#     # best_rads_dict["Version"].append("4_closest")
#     # best_rads_dict["Project"].append(project)
#     # best_rads_dict["Radius"].append(best_rad)
#
#     # best radius, 4 optimal
#     rad_means = temp_df.query("Condition=='optimal'").groupby("Radius")["Magnitude"].mean()
#     best_rad = rad_means[rad_means>0.3].idxmin()
#     mag_axes[proj_idx].plot(best_rad, rad_means[best_rad], marker='o',
#                             color="r", markersize=10)
#     best_rads_dict["Version"].append("4_optimal")
#     best_rads_dict["Project"].append(project)
#     best_rads_dict["Radius"].append(best_rad)
#
# # because P6 is missing
# mag_axes[-1].axis("off")
# foc_axes[-1].axis("off")
# mag_fig.savefig(join(fig_dir, "CvsO_mag_projavg.png"))
# foc_fig.savefig(join(fig_dir, "CvsO_foc_projavg.png"))
#
# rad_df = pd.DataFrame.from_dict(best_rads_dict)
#
# data_dir = join(root_dir, "simnibs/3_results")
# df_3 = pd.read_pickle(join(data_dir, "df_3.pickle"))
# df_3["Version"] = pd.Series(np.ones(len(df_3))*3)
# data_dir = join(root_dir, "simnibs/4_results")
# df_4 = pd.read_pickle(join(data_dir, "df_4.pickle"))
# df_4["Version"] = pd.Series(np.ones(len(df_4))*4)
# df = pd.concat([df_3, df_4])
# df_dict = {"Subject":[], "Radius":[], "Magnitude":[], "Focality":[],
#            "Project":[], "Version":[]}
# versions = ["3_closest", "4_closest", "4_optimal"]
# for proj_idx, project in enumerate(projects):
#     for version in versions:
#         vers = 4. if "4" in version else 3.
#         cond = version[2:]
#         this_df = df.query(f"Project=='{project}' and Version=={vers} and "
#                            f"Condition=='{cond}'")
#         this_radius = rad_df.query(f"Version=='{version}' and "
#                                    f"Project=='{project}'")["Radius"].values[0]
#
#         # TEMPORARY
#         this_radius = 50 if this_radius < 50 else this_radius
#
#         # build a new df in a format where this can be easily plotted
#         radii = this_df.iloc[0]["Radii"]
#         rad_idx = np.where(radii==this_radius)[0]
#         for row_idx, row in this_df.iterrows():
#             df_dict["Subject"].append(row["Subject"])
#             df_dict["Radius"].append(this_radius)
#             df_dict["Magnitude"].append(row["Mags"][rad_idx][0])
#             df_dict["Focality"].append(row["Focs"][rad_idx][0])
#             df_dict["Project"].append(project)
#             df_dict["Version"].append(version)
#
# temp_df = pd.DataFrame.from_dict(df_dict)
#
# # magnitude
#
# fig, ax = plt.subplots(1, figsize=(38.4, 8))
# sns.violinplot(data=temp_df, x="Project", y="Magnitude", hue="Version", ax=ax,
#                inner="points", hue_order=versions)
# x_labels = []
# for proj in projects:
#     rads = []
#     for version in versions:
#         rads.append(temp_df.query(f"Project=='{proj}' and Version=='{version}'")["Radius"].values[0])
#     x_labels.append(f"{proj} ({rads[0]}/{rads[1]}/{rads[2]})")
# ax.set_xticklabels(x_labels)
# ax.axhline(0.2, color="blue", linestyle='--')
# ax.axhline(0.3, color="black", linestyle='--')
# fig.savefig(join(fig_dir, "adapted_radius_mag.png"))
#
# # focality
#
# fig, ax = plt.subplots(1, figsize=(38.4, 8))
# sns.violinplot(data=temp_df, x="Project", y="Focality", hue="Version", ax=ax,
#                inner="points", hue_order=versions)
# x_labels = []
# for proj in projects:
#     rads = []
#     for version in versions:
#         rads.append(rad_df.query(f"Project=='{proj}' and Version=='{version}'")["Radius"].values[0])
#     x_labels.append(f"{proj} ({rads[0]}/{rads[1]}/{rads[2]})")
# ax.set_xticklabels(x_labels)
# fig.savefig(join(fig_dir, "adapted_radius_foc.png"))

# 3 only
best_rads_dict = {"Project":[], "Radius":[]}
data_dir = join(root_dir, "simnibs/3_results")
dff = pd.read_pickle(join(data_dir, "df_3.pickle"))
df = dff.query("Condition=='closest'")
projects = np.sort(df["Project"].unique())
mag_fig, mag_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
plt.suptitle("Magnitude")
foc_fig, foc_axes = plt.subplots(2, 4, figsize=(38.4, 21.6))
plt.suptitle("Focality")
mag_axes = [ax for axe in mag_axes for ax in axe]
foc_axes = [ax for axe in foc_axes for ax in axe]
for proj_idx, project in enumerate(projects):
    this_df = df.query(f"Project=='{project}'")
    # build a new df in a format where this can be easily plotted
    df_dict = {"Subject":[], "Radius":[], "Magnitude":[], "Focality":[]}
    radii = this_df.iloc[0]["Radii"]
    for row_idx, row in this_df.iterrows():
        for rad_idx, radius in enumerate(radii):
            df_dict["Subject"].append(row["Subject"])
            df_dict["Radius"].append(radius)
            df_dict["Magnitude"].append(row["Mags"][rad_idx])
            df_dict["Focality"].append(row["Focs"][rad_idx])
    temp_df = pd.DataFrame.from_dict(df_dict)

    sns.lineplot(data=temp_df, x="Radius", y="Magnitude",
                 ax=mag_axes[proj_idx])
    sns.lineplot(data=temp_df, x="Radius", y="Focality",
                 ax=foc_axes[proj_idx])
    mag_axes[proj_idx].set_title(project)
    mag_axes[proj_idx].set_ylim([.05, .7])
    foc_axes[proj_idx].set_ylim([1000, 9000])
    foc_axes[proj_idx].set_title(project)

    # best radius, 3
    rad_means = temp_df.groupby("Radius")["Magnitude"].mean()
    best_rad = rad_means[rad_means>0.2].idxmin()
    mag_axes[proj_idx].plot(best_rad, rad_means[best_rad], marker='o',
                            color="black", markersize=10)
    mag_axes[proj_idx].axhline(0.2, color="black", linestyle="--")
    best_rads_dict["Project"].append(project)
    best_rads_dict["Radius"].append(best_rad)

# because P6 is missing
mag_axes[-1].axis("off")
foc_axes[-1].axis("off")
mag_fig.savefig(join(fig_dir, "rads_3only_mag.pdf"))
foc_fig.savefig(join(fig_dir, "rads_3only_foc.pdf"))

rad_df = pd.DataFrame.from_dict(best_rads_dict)

df_dict = {"Subject":[], "Radius":[], "Magnitude":[], "Focality":[],
           "Project":[]}
for proj_idx, project in enumerate(projects):
    this_df = df.query(f"Project=='{project}'")
    this_radius = rad_df.query(f"Project=='{project}'")["Radius"].values[0]

    # build a new df in a format where this can be easily plotted
    radii = this_df.iloc[0]["Radii"]
    rad_idx = np.where(radii==this_radius)[0]
    for row_idx, row in this_df.iterrows():
        df_dict["Subject"].append(row["Subject"])
        df_dict["Radius"].append(this_radius)
        df_dict["Magnitude"].append(row["Mags"][rad_idx][0])
        df_dict["Focality"].append(row["Focs"][rad_idx][0])
        df_dict["Project"].append(project)

temp_df = pd.DataFrame.from_dict(df_dict)

# magnitude

fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=temp_df, x="Project", y="Magnitude", ax=ax,
               inner=None)
sns.stripplot(data=temp_df, x="Project", y="Magnitude", ax=ax,
              size=15, color="black")
x_labels = []
for proj in projects:
    rad = temp_df.query(f"Project=='{proj}'")["Radius"].values[0]
    x_labels.append(f"{proj} ({rad})")
ax.set_xticklabels(x_labels)
ax.axhline(0.2, color="black", linestyle='--')
fig.savefig(join(fig_dir, "adapted_radius_mag_3only.pdf"))

# focality

fig, ax = plt.subplots(1, figsize=(38.4, 8))
sns.violinplot(data=temp_df, x="Project", y="Focality", ax=ax, inner=None)
sns.stripplot(data=temp_df, x="Project", y="Focality", ax=ax,
              size=15, color="black")
x_labels = []
for proj in projects:
    rad = rad_df.query(f"Project=='{proj}'")["Radius"].values[0]
    x_labels.append(f"{proj} ({rad})")
ax.set_xticklabels(x_labels)
fig.savefig(join(fig_dir, "adapted_radius_foc_3only.pdf"))



# # plot by subject for a certain version
#
# version = 3
# data_dir = join(root_dir, f"simnibs/{version}_results")
# df = pd.read_pickle(join(data_dir, f"df_{version}.pickle"))
#
# projects = np.sort(df["Project"].unique())
#
# for project in projects:
#     this_df = df.query(f"Project=='{project}'")
#     subjs = np.sort(this_df["Subject"].unique())
#
#     fig = plt.figure(figsize=(38.4, 21.6))
#     plt.suptitle(f"{project} (v{version})", fontsize=28)
#     col_n = 10
#     if version == 4:
#         row_n = 4
#         valid_rows = 2
#         alpha = 0.5
#     else:
#         row_n = 5
#         valid_rows = 1
#         alpha = 1.
#
#     outer = gridspec.GridSpec(row_n, col_n)
#     for subj_idx, subj in enumerate(subjs):
#         subj_df = this_df.query(f"Subject=='{subj}'")
#         if len(subj_df) != valid_rows:
#             continue
#         inner = gridspec.GridSpecFromSubplotSpec(2, 1,
#                                                  subplot_spec=outer[subj_idx])
#
#         closest_radii = subj_df.query("Condition=='closest'")["Radii"].values[0]
#         closest_mags = subj_df.query("Condition=='closest'")["Mags"].values[0]
#         closest_focs = subj_df.query("Condition=='closest'")["Focs"].values[0]
#         if version == 4:
#             optimal_radii = subj_df.query("Condition=='optimal'")["Radii"].values[0]
#             optimal_mags = subj_df.query("Condition=='optimal'")["Mags"].values[0]
#             optimal_focs = subj_df.query("Condition=='optimal'")["Focs"].values[0]
#
#             if not np.array_equal(closest_radii, optimal_radii):
#                 raise ValueError("Radii of closest and optimal not the same.")
#             else:
#                 radii = closest_radii
#         else:
#             radii = closest_radii
#
#         ax = plt.Subplot(fig, inner[0])
#         ax.set_title(subj)
#         ax.plot(radii, closest_mags, label="closest", alpha=alpha)
#         if version == 4:
#             ax.plot(radii, optimal_mags, label="optimal", alpha=alpha)
#         ax.set_ylim([0.05, 1.])
#         ax.set_xticks(radii)
#         ax.set_xticklabels([])
#         ax.axhline(0.2, linestyle="--", color="black", alpha=.5)
#         if outer[subj_idx].colspan[0] == 0:
#             ax.set_ylabel("Mag.", fontsize=18, fontweight="bold")
#             ax.set_yticks(np.round(np.arange(0.1, 1., .2), 1))
#         else:
#             ax.set_yticks([])
#         fig.add_subplot(ax)
#
#         ax = plt.Subplot(fig, inner[1])
#         ax.plot(radii, closest_focs, label="closest", alpha=.5)
#         if version == 4:
#             ax.plot(radii, optimal_focs, label="optimal", alpha=.5)
#         ax.set_ylim([250, 5000])
#         if outer[subj_idx].colspan[0] == 0:
#             ax.set_ylabel("Foc.", fontsize=18, fontweight="bold")
#             ax.set_yticks(np.arange(500, 7500, 2000))
#         else:
#             ax.set_yticks([])
#         ax.set_xticks(radii)
#         if outer[subj_idx].rowspan[0] != (row_n - 1):
#             ax.set_xticklabels([])
#         else:
#             ax.set_xlabel("Radii", fontsize=18, fontweight="bold")
#         fig.add_subplot(ax)
#     handles, labels = ax.get_legend_handles_labels()
#     fig.legend(handles, labels, loc="upper right")
#     plt.savefig(join(fig_dir, f"focmag_{project}_v{version}.png"))
