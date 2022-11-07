import numpy as np
import os
import shutil
from scipy.io import savemat
from time import perf_counter
from datetime import timedelta

from simnibs import __version__, sim_struct, mesh_io, mni2subject_coords
from simnibs.utils.file_finder import SubjectFiles
import Nx1_stuff
from emp_chandefs import prepare_emp

def emp_montage(subj_dict, proj_dict, root_dir, extract_only=False):
    """Simulate previous literature montages"""
    project, mask, hemi = (list(proj_dict.keys())[0],
                           *list(proj_dict.values())[0])
    print(f"\n\n\n{project} {mask} {hemi}\n\n\n")
    begin_time = perf_counter()
    subname, subpath = list(subj_dict.keys())[0], list(subj_dict.values())[0]
    version = int(__version__[0])
    if version > 3:
        var_name = 'E_magn'
        field_name = "magnE"
    else:
        var_name = 'E_norm'
        field_name = "normE"
    subject_files = SubjectFiles(subpath=subpath)
    pathfem = os.path.join(root_dir, f"{version}_emp",
                           f"{subname}_{project}")
    if os.path.isdir(pathfem) and not extract_only:
        print("Already exists. Skipping.")
        return None
    print(pathfem)
    print(f"Subject {subject_files.subid}")

    try:
        m = mesh_io.read_msh(subject_files.fnamehead)
    except:
        print("No mesh file found.")
        return None
    if version > 3:
        m = Nx1_stuff.relabel_internal_air(m, subpath)

    if mask:
        mask_path = os.path.join(root_dir, "ROI", mask)
        _, mask_pos = Nx1_stuff.convert_mask(mask_path, hemi, subpath)
        pos_center = Nx1_stuff.get_closest_skin_pos(mask_pos, m)

        pathfem = os.path.abspath(os.path.expanduser(pathfem))
        if not os.path.isdir(pathfem):
            os.mkdir(pathfem)
        if not extract_only:
            mesh_io.write_geo_spheres([pos_center],
                                       os.path.join(pathfem,
                                                    'mesh_with_ROI.geo'),
                                       name=('center'))
            mesh_io.write_msh(m, os.path.join(pathfem, 'mesh_with_ROI.msh'))

    if not extract_only:
        S = prepare_emp(project)
        S.subpath = subpath
        if project == "P2" or project == "P2_5050" or project == "P6":
            S.eeg_cap = S.subpath + '/eeg_positions' + '/EEGcap_incl_cheek_buci_3.csv'
        S.pathfem = pathfem
        S.map_to_surf = True
        S.map_to_vol = True
        S.map_to_fsavg = True
        S.map_to_MNI = True
        S.open_in_gmsh = False
        S.fnamehead = subject_files.fnamehead
        S.run()

    if project == "P6":
        # P6 has a spherical ROI, not cortical mask
        msh_file = "TDCS_1_scalar.msh"
        msh_file = "T1w.nii_" + msh_file if version > 3 else subject_files.subid + "_" + msh_file

        mesh = mesh_io.read_msh(os.path.join(pathfem, msh_file))
        gray_matter = mesh.crop_mesh(2)
        ROI_center = [13, -79, -37]
        subj_center = mni2subject_coords(ROI_center, subpath)
        rad = 10.
        elm_centers = gray_matter.elements_baricenters()[:]
        roi = np.linalg.norm(elm_centers - subj_center, axis=1) < rad
        elm_vols = gray_matter.elements_volumes_and_areas()[:]
        gray_matter.add_element_field(roi, 'roi')
        field = gray_matter.field[field_name][:]
        vals = field[roi]
        mean = np.average(vals, weights=elm_vols[roi])
        median = np.median(vals)
        focality_med = np.sum(vals[vals>median])
        focality_mean = np.sum(vals[vals>mean])
        if not extract_only:
            mesh_io.write_msh(gray_matter, os.path.join(S.pathfem,
                                                        "results.msh"))
        pos_center = subj_center
    else:
        msh_file = "TDCS_1_scalar_central.msh"
        msh_file = "T1w.nii_" + msh_file if version > 3 else subject_files.subid + "_" + msh_file
        m_surf = Nx1_stuff.get_central_gm_with_mask(subpath, hemi, mask_path)
        nd_sze = m_surf.nodes_volumes_or_areas().value
        idx_mask = m_surf.nodedata[0].value
        try:
            m = mesh_io.read_msh(os.path.join(pathfem, "subject_overlays",
                                              msh_file))
        except:
            return None
        assert m.nodes.nr == m_surf.nodes.nr
        nd = next(x.value for x in m.nodedata if x.field_name==var_name)
        m_surf.add_node_field(nd, "result")
        median = np.median(nd[idx_mask])
        mean = np.mean(nd[idx_mask])
        focality_med = np.sum(nd_sze[nd > median])
        focality_mean = np.sum(nd_sze[nd > mean])
        if not extract_only:
            mesh_io.write_msh(m_surf, os.path.join(pathfem, 'results.msh'))

    mdic = {"pos_center": pos_center,
            "focality_med": focality_med,
            "focality_mean": focality_mean,
            "median": median,
            "mean": mean
            }
    savemat(os.path.join(pathfem, 'summary_metrics.mat'), mdic)

def rad_only(subj_dict, mask_dict, condition, radii, EL_center,
             EL_surround, root_dir, N=3, cutoff=.1,
             multichannel=True, current_center=0.002):
    """Simulate with variable radius montages"""
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
            "best_phi": phi,
            "roi_median_f": roi_median_f,
            "focality_f": focality_f
            }
    savemat(os.path.join(pathfem, 'summary_metrics.mat'), mdic)
    end_time = perf_counter()
    time_str = str(timedelta(seconds=(end_time - begin_time)))
    print(f"\n{subject_files.subid} finished in {time_str}.\n")
