from utils import *
import matplotlib.pyplot as plt
import os
import re

# detect whether running in HGW or Berlin
if os.path.isdir("/home/hannaj"):
    root_dir = "/media/Linux1_Data02/antonenkod/AD-Stim_SimNibs4/"
    fig_dir = f"/home/hannaj/temp/emp_temp/"
else:
    root_dir = "/home/jev/simnibs/"
    fig_dir = f"/home/jev/temp/emp_temp/"


overwrite = False
# get subject names
subjs = os.listdir(f"{root_dir}")

cam_dist = 400
fig_count = 0
fig_sub_idx = 0
for subj_idx, subj in enumerate(subjs):
    if "TS3-MCI" not in subj:
        continue
    if fig_sub_idx == 0:
        fig, axes = plt.subplots(3, 8, figsize=(38.4, 12.))
    # load mesh
    mesh = pv.read(f"{root_dir}{subj}/simu_F3-AF4_5x7/{subj}_TDCS_1_scalar.msh")
    # plot fields
    mag_image = mag_plot(mesh, cam_dist=cam_dist, foc="elec_emp")
    axes[2, fig_sub_idx].imshow(mag_image)
    axes[2, fig_sub_idx].axis("off")

    # plot electrodes
    elec_image, foc = elec_plot_emp(mesh, return_foc=True,
                                    cam_dist=cam_dist)
    axes[0, fig_sub_idx].imshow(elec_image)
    axes[0, fig_sub_idx].axis("off")

    elec_image, foc = elec_plot_emp(mesh, return_foc=True,
                                    cam_dist=cam_dist, foc_elec="cathode")
    axes[1, fig_sub_idx].imshow(elec_image)
    axes[1, fig_sub_idx].axis("off")

    # text
    trans = axes[0, fig_sub_idx].transAxes
    axes[0, fig_sub_idx].text(0.1, 0.8, subj[-3:], fontsize=36,
                           transform=trans)
    fig_sub_idx += 1

    if fig_sub_idx == 8 or (subj_idx==(len(subjs)-1)):
        plt.tight_layout()
        plt.savefig(f"{fig_dir}{fig_count}_magn.pdf")
        plt.close()
        fig_count += 1
        fig_sub_idx = 0
