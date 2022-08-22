import numpy as np
import os
import shutil
from scipy.io import savemat
from time import perf_counter
from datetime import timedelta

from simnibs import __version__, sim_struct, mesh_io
from simnibs.utils.file_finder import SubjectFiles
import Nx1_stuff

def build_subject_paths(root_dir, version):
    subjectpaths = []
    subj_dirs = next(os.walk(root_dir))[1]
    for subj_dir in subj_dirs:
        if round(version) == 3:
            m2m_str = f"m2m_T1w"
        elif round(version) == 4:
            m2m_str = f"m2m_T1w.nii"
        subjectpaths.append(os.path.join(root_dir, subj_dir, m2m_str))
    return subjectpaths, subj_dirs

mask = 'P1_rTP-RH'
phi = 35.
hemi = 'rh'
condition = 'closest' # 'closest' or 'optimal'
version = int(__version__[0])

root_dir = "/media/Linux5_Data03/hannaj/simnibs/"
data_dir = os.path.join(root_dir, str(round(version)))
subjectpaths, subj_dirs = build_subject_paths(data_dir, version)

print(f"\nVersion {version}\n")

if int(version)>3:
    var_name = 'E_magn'
else:
    var_name = 'E_norm'
cutoff = 0.2 # in [V/m]

N = 3 # number of surround electrodes
radius_surround = [60., 65., 70., 75., 80., 85., 90., 95.] # in [mm]
multichannel = True # set to True if current in each surround electrode is separately controlled
current_center = 0.002  # Current flow through center channel (A)

# properties of centre electrode
EL_center = sim_struct.ELECTRODE()
EL_center.shape = 'ellipse'  # round shape
EL_center.dimensions = [20, 20]  # 30 mm diameter
EL_center.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer

# properties of surround electrodes
EL_surround = sim_struct.ELECTRODE()
EL_surround.shape = 'ellipse'  # round shape
EL_surround.dimensions = [20, 20]  # 30 mm diameter
EL_surround.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer


###  PREPARE SIMULATIONS
############################
assert condition == 'closest' or condition == 'optimal'
for subpath, subj, in zip(subjectpaths, subj_dirs):
    begin_time = perf_counter()
    subject_files = SubjectFiles(subpath=subpath)
    pathfem = os.path.join(root_dir, f"{version}_results",
                           f"{mask}__{subj}__{condition}")
    print(pathfem)
    print(f"Subject {subject_files.subid}")

    # load head mesh
    try:
        m = mesh_io.read_msh(subject_files.fnamehead)
    except:
        print("No mesh found. Skipping...")
        continue
    if int(__version__[0])>3:
        m = Nx1_stuff.relabel_internal_air(m, subpath)

    # convert mask to individual space (on central surface)
    mask_path = os.path.join(root_dir, "ROI", mask)
    _, mask_pos = Nx1_stuff.convert_mask(mask_path, hemi, subpath)
    if condition == 'closest':
        # use skin position closest to CoG of mask
        pos_center = Nx1_stuff.get_closest_skin_pos(mask_pos, m)
    else:
        # project mask positions to pial surface of tet mesh
        # and relabel corresponding GM triangles
        m = Nx1_stuff.project_to_pial(mask_pos, m)
        # solve FEM to get optimal position on skin
        # with lowest ohmic ressitance to mask
        m, pos_center = Nx1_stuff.get_optimal_center_pos(m)
    EL_center.centre = pos_center

    # write out mesh with ROI and position of central electrode as control
    pathfem = os.path.abspath(os.path.expanduser(pathfem))
    if not os.path.isdir(pathfem):
        os.mkdir(pathfem)
    mesh_io.write_geo_spheres([pos_center],
                              os.path.join(pathfem,'mesh_with_ROI.geo'),
                                           name=('center'))
    mesh_io.write_msh(m, os.path.join(pathfem,'mesh_with_ROI.msh'))


    ###  RUN SIMULATIONS FOR VARIING RADII
    #######################################
    try:
        Nx1_stuff.run_simus(subpath, os.path.join(pathfem,'radius'),
                        current_center, N, radius_surround, [0.],
                        EL_center, EL_surround)
    except:
        print("Skipping...")
        continue

    try:
        m_surf, roi_median_r, focality_r, best_radius = Nx1_stuff.analyse_simus(subpath,
                                            os.path.join(pathfem,'radius'),
                                            hemi, mask_path,
                                            radius_surround, [phi],
                                            var_name, cutoff)
    except:
        print("Skipping...")
        continue

    best_radius = best_radius[0]
    mesh_io.write_msh(m_surf,os.path.join(pathfem,'results_radii.msh'))
    print(radius_surround)
    print(roi_median_r)
    print(focality_r)
    print('selecting radius '+ str(best_radius))


    ### RUN FINAL SIMULATION AND SAFE ...
    ###############################################
    Nx1_stuff.run_simus(subpath, os.path.join(pathfem,'final'),
                        current_center, N, [best_radius], [phi],
                        EL_center, EL_surround)

    m_surf, roi_median_f, focality_f, _ = Nx1_stuff.analyse_simus(subpath,
                                            os.path.join(pathfem,'final'),
                                            hemi, mask_path,
                                            [best_radius], [phi],
                                            var_name, cutoff)
    mesh_io.write_msh(m_surf,os.path.join(pathfem,'results_final.msh'))

    mdic = {"pos_center": pos_center,
            "radius_surround": radius_surround,
            "roi_median_r": roi_median_r,
            "focality_r": focality_r,
            "best_radius": best_radius,
            "phi_offset": 0,
            "roi_median_p": 0,
            "focality_p": 0,
            "best_phi": phi,
            "final_radius": best_radius,
            "roi_median_f": roi_median_f,
            "focality_f": focality_f
            }
    savemat(os.path.join(pathfem, 'summary_metrics.mat'), mdic)
    end_time = perf_counter()
    time_str = str(timedelta(seconds=(end_time - begin_time)))
    print(f"\n{subject_files.subid} finished in {time_str}.\n")
