from utils import *
import matplotlib.pyplot as plt
import os
import re

projects = [f"P{idx}" for idx in range(1,9)]

# detect whether running in HGW or Berlin
if os.path.isdir("/home/hannaj"):
    root_dir = "/home/hannaj/simnibs/"
else:
    root_dir = "/home/jev/simnibs/"
res_dir = f"{root_dir}4_emp/"
fig_dir = f"{res_dir}figures/"

overwrite = False


# get subject names
subjs = os.listdir(f"{root_dir}4/")

cam_dist = 400

for subj in subjs:
    outfile = f"{fig_dir}{subj}_magn.pdf"
    if f"{subj}_magn.pdf" in os.listdir(fig_dir) and not overwrite:
        print(f"{outfile} already exists. Skipping...")
        continue
    fig, axes = plt.subplots(3, 8, figsize=(38.4, 12.))
    for proj_idx, project in enumerate(projects):
        # load mesh
        mesh = pv.read(f"{res_dir}{subj}_{project}/{subj}_TDCS_1_scalar.msh")
        # plot fields
        mag_image = mag_plot(mesh, cam_dist=cam_dist, foc="elec_emp")
        axes[2, proj_idx].imshow(mag_image)
        axes[2, proj_idx].axis("off")

        # plot electrodes
        P8 = False
        if project == "P8":
            P8 = True
        elec_image, foc = elec_plot_emp(mesh, return_foc=True,
                                        cam_dist=cam_dist, P8=P8)
        axes[0, proj_idx].imshow(elec_image)
        axes[0, proj_idx].axis("off")
        if not (project == "P7" or project == "P8"):
            elec_image, foc = elec_plot_emp(mesh, return_foc=True,
                                            cam_dist=cam_dist, P8=P8,
                                            foc_elec="cathode")
        axes[1, proj_idx].imshow(elec_image)
        axes[1, proj_idx].axis("off")

        # text
        trans = axes[0, proj_idx].transAxes
        axes[0, proj_idx].text(0.8, 0.1, project, fontsize=36,
                               transform=trans)

        trans = axes[0, 0].transAxes
        axes[0, 0].text(0.005, 0.8, subj, fontsize=58, transform=trans)

    plt.tight_layout()
    plt.savefig(f"{fig_dir}{subj}_magn.pdf")
    plt.close()
