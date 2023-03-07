from utils import *
import matplotlib.pyplot as plt
import os

masks = ["P1_rTP-RH", "P2_lPCC-new-LH",  "P3_lTP-LH", "P4_lIFG-LH", "P5_lM1-LH",
         "P7_rDLPFCnew-RH", "P8_lDLPFC-LH"]
if os.path.isdir("/home/hannaj"):
    root_dir = "/home/hannaj/simnibs/"
else:
    root_dir = "/home/jev/simnibs/"
res_dir = f"{root_dir}4_results/"
fig_dir = f"{root_dir}figures/"
algos = ["closest"]

# get subject names
subjs = os.listdir(f"{root_dir}4/")

cam_dist = 400

for subj in subjs:
    #try:
    fig, axes = plt.subplots(len(algos)*2+1, len(masks),
                             figsize=(38.4, 12.))
    for mask_idx, mask in enumerate(masks):
        for alg_idx, alg in enumerate(algos):
            subj_dir = f"{res_dir}{mask}__{subj}__{alg}/final/"
            # load mesh
            mesh = pv.read(f"{subj_dir}{subj}_TDCS_1_scalar.msh")

            mag_image = mag_plot(mesh, cam_dist=cam_dist)
            axes[alg_idx*len(algos)+2, mask_idx].imshow(mag_image)
            axes[alg_idx*len(algos)+2, mask_idx].axis("off")

            elec_image, foc = elec_plot(mesh, return_foc=True,
                                        cam_dist=cam_dist)
            axes[alg_idx*2+1, mask_idx].imshow(elec_image)
            axes[alg_idx*2+1, mask_idx].axis("off")

        roi_dir = f"{res_dir}{mask}__{subj}__closest/"
        roi_mesh = pv.read(f"{roi_dir}results_final.msh")
        roi_image = roi_plot(roi_mesh, foc=foc, cam_dist=cam_dist)
        axes[0, mask_idx].imshow(roi_image)
        axes[0, mask_idx].axis("off")

        # text
        trans = axes[0, mask_idx].transAxes
        axes[0, mask_idx].text(0.8, 0.1, mask[:2], fontsize=36,
                               transform=trans)

    trans = axes[1, 0].transAxes
    axes[1, 0].text(0.01, 0.9, "Closest", fontsize=24, transform=trans)
    # trans = axes[3, 0].transAxes
    # axes[3, 0].text(0.01, 0.9, "Optimal", fontsize=24, transform=trans)

    trans = axes[0, 0].transAxes
    axes[0, 0].text(0.005, 0.8, subj, fontsize=58, transform=trans)

    plt.tight_layout()
    plt.savefig(f"{fig_dir}{subj}_magn.pdf")
    plt.close()
    # except:
    #     print(f"Could not process {subj}.")
