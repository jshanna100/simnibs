import numpy as np
import os
import shutil
from scipy.io import savemat
from time import perf_counter
from datetime import timedelta

from simnibs import __version__, sim_struct, mesh_io
from simnibs.utils.file_finder import SubjectFiles
import Nx1_stuff


def rad_only(subj_dict, mask_dict, condition, radii, EL_center,
             EL_surround, root_dir, N=3, cutoff=.1,
             multichannel=True, current_center=0.002):
    begin_time = perf_counter()
    subname, subpath = list(subj_dict.keys())[0], list(subj_dict.values())[0]
    mask, phi, hemi = (list(mask_dict.keys())[0], *list(mask_dict.values())[0])
    version = int(__version__[0])
    if int(version)>3:
        var_name = 'E_magn'
    else:
        var_name = 'E_norm'

    subject_files = SubjectFiles(subpath=subpath)
    pathfem = os.path.join(root_dir, f"{version}_results",
                           f"{mask}__{subname}__{condition}")
    if os.path.isdir(pathfem):
        print("Already exists. Skipping.")
        return None
    print(pathfem)
    print(f"Subject {subject_files.subid}")

    try:
        m = mesh_io.read_msh(subject_files.fnamehead)
    except:
        print("No mesh file found.")
        return None
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
                        current_center, N, radii, [0.],
                        EL_center, EL_surround)

        m_surf, roi_median_r, focality_r, best_radius = Nx1_stuff.analyse_simus(subpath,
                                            os.path.join(pathfem,'radius'),
                                            hemi, mask_path,
                                            radii, [phi],
                                            var_name, cutoff)
    except:
        print("Could not process.")
        return None

    best_radius = best_radius[0]
    mesh_io.write_msh(m_surf,os.path.join(pathfem,'results_radii.msh'))
    print(radii)
    print(roi_median_r)
    print(focality_r)
    print('selecting radius '+ str(best_radius))


    ### RUN FINAL SIMULATION AND SAFE ...
    ###############################################
    try:
        Nx1_stuff.run_simus(subpath, os.path.join(pathfem,'final'),
                            current_center, N, [best_radius], [phi],
                            EL_center, EL_surround)

        m_surf, roi_median_f, focality_f, _ = Nx1_stuff.analyse_simus(subpath,
                                                os.path.join(pathfem,'final'),
                                                hemi, mask_path,
                                                [best_radius], [phi],
                                                var_name, cutoff)
        mesh_io.write_msh(m_surf,os.path.join(pathfem,'results_final.msh'))
    except:
        print("Could not process.")
        return None

    mdic = {"pos_center": pos_center,
            "radius_surround": radii,
            "roi_median_r": roi_median_r,
            "focality_r": focality_r,
            "best_radius": best_radius,
            "phi_offset": 0, # no longer relevant with constant phi
            "roi_median_p": 0,  # no longer relevant with constant phi
            "focality_p": 0,  # no longer relevant with constant phi
            "best_phi": phi,
            "final_radius": best_radius,  # same as best with constant phi
            "roi_median_f": roi_median_f,
            "focality_f": focality_f
            }
    savemat(os.path.join(pathfem, 'summary_metrics.mat'), mdic)
    end_time = perf_counter()
    time_str = str(timedelta(seconds=(end_time - begin_time)))
    print(f"\n{subject_files.subid} finished in {time_str}.\n")
