from project_settings import projects
from utils import (elec_plot_geo, pvmesh_from_skin_geo, elpos_from_geo,
                   elpos_from_brainsight)
from os.path import join
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

proj_nr = 1
subj = "ernie"
cam_dist = 400
base_dir = "/home/jev/simnibs/"

proj = projects[proj_nr]
radius = proj.radius[0]
phi = proj.phi

mesh_dir = join(base_dir, "test_closest", f"P{proj_nr}_{subj}")
mesh = pvmesh_from_skin_geo(join(mesh_dir, f"P{proj_nr}_{subj}_skin.geo"))
points = elpos_from_geo(join(mesh_dir,
                             f"P{proj_nr}_{subj}_{int(radius)}_elpos.geo"))
elec_img = elec_plot_geo(mesh, points, return_foc=False, cam_dist=cam_dist)
ras_coords = elpos_from_brainsight(join(mesh_dir,
                                        f"P{proj_nr}_{subj}_{int(radius)}brainsight_RAS.txt"))
lps_coords = elpos_from_brainsight(join(mesh_dir,
                                        f"P{proj_nr}_{subj}_{int(radius)}brainsight_LPS.txt"))

colors = ["black", "red", "green", "blue"]
y_base = 0.99
x_base = 0.01
fig, axes = plt.subplots(2, 1, figsize=(15.32, 21.6))
kwargs = {"transform":axes[0].transAxes, "verticalalignment":"top"}
axes[0].text(x_base, y_base, "MeMoSLAP Placement Guide\n",
             **kwargs, weight="bold", fontsize=38)
axes[0].text(x_base, y_base-0.1,
             f"Project P{proj_nr}\nID: {subj}\nRadius: {int(radius)}mm",
             **kwargs, fontsize=28)
axes[0].text(x_base, y_base-0.35, "Coordinates (XYZ)", **kwargs, fontsize=28, weight="bold")
axes[0].text(x_base, y_base-0.45, "RAS", **kwargs, fontsize=28)
for idx, ((k, v), color) in enumerate(zip(ras_coords.items(), colors)):
    k = k.replace("_electrode", "")
    k = k.replace("surround", "rad")
    v = np.round(v, 2)
    axes[0].text(x_base, y_base-0.55-idx*0.1, f"{k}: {v[0]:.2f}, {v[1]}, {v[2]}",
                 **kwargs, fontsize=28, color=color)

x_base = 0.51
axes[0].text(x_base, y_base-0.45, "LPS", **kwargs, fontsize=28)
for idx, ((k, v), color) in enumerate(zip(lps_coords.items(), colors)):
    k = k.replace("_electrode", "")
    k = k.replace("surround", "rad")
    v = np.round(v, 2)
    axes[0].text(x_base, y_base-0.55-idx*0.1, f"{k}: {v[0]:.2f}, {v[1]}, {v[2]}",
                 **kwargs, fontsize=28, color=color)
axes[0].axis("off")

axes[1].imshow(elec_img)
axes[1].axis("off")

plt.tight_layout()

plt.savefig(join(mesh_dir, f"placeguide_P{proj_nr}_{subj}.pdf"))
