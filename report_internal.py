from project_settings import projects
from utils import (elec_plot_geo, pvmesh_from_skin_geo, elpos_from_geo,
                   elpos_from_brainsight, mag_plot, roi_plot)
from os.path import join
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import pickle
plt.ion()

proj_nr = 1
subj = "ernie"
cam_dist = 400
base_dir = "/home/jev/simnibs/"

proj = projects[proj_nr]
radius = proj.radius[0]
phi = proj.phi

mesh_dir = join(base_dir, "test_closest", f"P{proj_nr}_{subj}")
with open(join(mesh_dir, "simnibs_memoslap_results.pkl"), "rb") as f:
    analysis = pickle.load(f)[5]
rad_idx = analysis["radius"].index(int(radius))
fieldmed = analysis["roi_median"]["E_magn"][rad_idx]
fieldsumsq = analysis["roi_squared"]["E_magn"][rad_idx]
foc = analysis["focality"]["E_magn"][rad_idx]

y_base = 0.99
x_base = 0.01
mos_str = """
          AB
          CD
          """
fig, axes = plt.subplot_mosaic(mos_str, figsize=(15.32, 21.6))
kwargs = {"transform":axes["A"].transAxes, "verticalalignment":"top"}
axes["A"].text(x_base, y_base, "MeMoSLAP FEM Summary",
             **kwargs, weight="bold", fontsize=38)
axes["A"].text(x_base, y_base-0.1,
             f"Project P{proj_nr}\n"
             f"ID: {subj}\n"
             f"Radius: {int(radius)}mm\n"
             f"Phi: {int(phi)}\n\n"
             f"Group Allocation: Target\n"
             f"Condition order: Target-Sham\n\n"
             f"Median field magnitude: {np.round(fieldmed, 3)}\n"
             f"Field magnitude sum²: {int(np.round(fieldsumsq))}\n"
             f"Focality: {int(np.round(foc))}\n\n"
             f"Baseline scan date: xx-xx-xxx\n"
             f"Simulated with version 4\n"
             f"  on xx-xx-xxxx\n"
             f"  by Jevri Hanna",
             **kwargs, fontsize=28)

axes["A"].axis("off")

# bottom panel
mesh = pvmesh_from_skin_geo(join(mesh_dir, f"P{proj_nr}_{subj}_skin.geo"))
points = elpos_from_geo(join(mesh_dir,
                             f"P{proj_nr}_{subj}_{int(radius)}_elpos.geo"))
elec_img, foc = elec_plot_geo(mesh, points, return_foc=True, cam_dist=cam_dist)
axes["B"].imshow(elec_img)
axes["B"].set_title("Electrodes", fontsize=30, weight="bold")
axes["B"].axis("off")

mesh = pv.read(join(mesh_dir, f"P{proj_nr}_{subj}_{int(radius)}.msh"))
mag_img = mag_plot(mesh, foc=foc, cam_dist=300)
axes["C"].imshow(mag_img)
axes["C"].set_title("Field Magnitude", fontsize=30, weight="bold")
axes["C"].axis("off")

roi_img = roi_plot(mesh, foc=foc, cam_dist=300)
axes["D"].imshow(roi_img)
axes["D"].set_title("ROI", fontsize=30, weight="bold")
axes["D"].axis("off")

#plt.tight_layout()
plt.savefig(join(mesh_dir, f"FEMSummary_P{proj_nr}_{subj}.pdf"))
