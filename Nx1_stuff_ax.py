# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 16:50:29 2021

@author: axthi
"""

import copy
import nibabel as nib
import numpy as np
import scipy.ndimage.morphology as mrph
import os
import re
from glob import glob, escape

from simnibs import mesh_io, sim_struct, __version__
from simnibs.simulation import cond, fem
from simnibs.utils.file_finder import SubjectFiles

if int(__version__[0])>3:
    from simnibs.utils import transformations
    from simnibs.utils.file_finder import get_reference_surf
else:
    from simnibs.msh import transformations
    from simnibs.utils.file_finder import templates


####################################
# start extra stuff for version 3
from copy import deepcopy
from scipy.spatial.transform import Rotation as R
from simnibs import file_finder

def sphereFit(pts, bounds = None):
    """
    Fit a circle or sphere to a point cloud.

    returns the radius and center points of the best fit sphere
    (adapted from https://jekel.me/2015/Least-Squares-Sphere-Fit/)

    Parameters
    ----------
    pts : array (Nx2) or (Nx3)
        point cloud

    Returns
    -------
    R : float64
        radius
    centre: ndarray (2,) or (3,)
        centre position

    """
    A = np.hstack((2*pts, np.ones((pts.shape[0],1)) ))
    f = np.sum(pts**2,1)
    C, residuals, rank, singval = np.linalg.lstsq(A,f,rcond=None)

    dim = pts.shape[1]
    R = np.sqrt(np.sum(C[0:dim]**2) + C[dim])
    return R, C[0:dim]


def sph2cart(az, el, r): # phi, theta, radius
    """Conversion from spherical to cartesian coordinates."""
    rcos_theta = r * np.cos(el)
    pts = np.zeros(( (3,) + rcos_theta.shape ))
    pts[0,:] = rcos_theta * np.cos(az)
    pts[1,:] = rcos_theta * np.sin(az)
    pts[2,:] = r * np.sin(el)
    return pts


def show_for_debugging(m,sph_centre,radius,P_centre,P_surround,surround_fit,M_sph):
    """Show some results in gmsh for debugging."""
    import tempfile
    fn_geo=tempfile.NamedTemporaryFile(suffix='.geo').name
    mesh_io.write_geo_spheres(sph_centre.reshape((1,3)),
                              fn_geo, name = 'center', mode = 'bw')
    mesh_io.write_geo_spheres(( M_sph @ [radius,0,0,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'x', mode = 'ba')
    mesh_io.write_geo_spheres(( M_sph @ [0,radius,0,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'y', mode = 'ba')
    mesh_io.write_geo_spheres(( M_sph @ [0,0,radius,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'z', mode = 'ba')

    TH_DEBUG = np.arange(-1.0,1.01,0.1)*np.pi
    PHI_DEBUG = np.arange(0.,1.01,0.05)*2*np.pi
    TH_DEBUG, PHI_DEBUG = np.meshgrid(TH_DEBUG, PHI_DEBUG)
    TH_DEBUG = TH_DEBUG.flatten()
    PHI_DEBUG = PHI_DEBUG.flatten()
    R_DEBUG = radius*np.ones_like(TH_DEBUG)

    pts=sph2cart(PHI_DEBUG, TH_DEBUG, R_DEBUG)
    pts = np.vstack(( pts, np.ones((1,pts.shape[1])) ))
    mesh_io.write_geo_spheres((M_sph @ pts)[:3,:].T, fn_geo,
                              name = 'sphere', mode = 'ba')

    mesh_io.write_geo_spheres( P_centre.reshape((1,3)),
                               fn_geo, name = 'centre', mode = 'ba')
    for i in range(len(P_surround)):
        mesh_io.write_geo_spheres( P_surround[i].reshape((1,3)),
                                  fn_geo, name = 'surr '+str(i), mode = 'ba')

    N_pts = 50
    for i in range(len(surround_fit)):
        tmp_centre = surround_fit[i][0]
        tmp_r = surround_fit[i][1]
        tmp_theta = surround_fit[i][2]
        tmp_theta_z0 = surround_fit[i][3]
        tmp_M = surround_fit[i][4]

        tmp_arc = np.vstack((
            tmp_r*np.sin(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)*np.arange(N_pts)/(N_pts-1)) + tmp_centre[0],
            np.zeros((1,N_pts)),
            tmp_r*np.cos(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)*np.arange(N_pts)/(N_pts-1)) + tmp_centre[1],
            np.ones((1,N_pts))
            ))
        tmp_arc=(tmp_M @ tmp_arc)[:3].T
        mesh_io.write_geo_spheres(tmp_arc, fn_geo, name = 'arc '+str(i), mode = 'ba')

    vis = mesh_io.gmsh_view.Visualization(m)
    vis.add_merge(fn_geo)
    vis.show()
    os.remove(fn_geo)


def get_surround_pos(center_pos, fnamehead, radius_surround=50, N=4,
                     pos_dir_1stsurround=None, phis_surround=None,
                     tissue_idx=1005, DEBUG=False):
    """
    Determine the positions of surround electrodes.

    Parameters
    ----------
    center_pos : array (3,) or string
        Center position of the central electrode.
    fnamehead : string
        Filename of head mesh.
    radius_surround : float or array (N,), optional
        Distance (centre-to-centre) between the centre and surround
        electrodes. The default is 50. Either a single number (same
        radius for all electrodes) or an array with N numbers
        (N: number of surround electrodes)
    N : int, optional
        Number of surround electrodes. The default is 4.
    pos_dir_1stsurround : array (3,) or string, optional
        A position indicating the direction from center_pos to
        the position of the first surround electrode. The default is None.
    phis_surround : array (N,), optional
        Angles in degree at which the electrodes will be place relative to the
        direction defined by pos_dir_1stsurround. The default is None, in which
        case the electrodes will be placed at [0, 1/N*360, ..., (N-1)/N*360]
        degrees.
    tissue_idx : int, optional
        Index of the tissue on which the surround positions will be planned on.
        (standard: 1005 for skin)
    DEBUG : Boolean, optional
        When set to True, a visualization in gmsh will open for control
        (standard: False)

    Returns
    -------
    P_surround : list of arrays (3,)
        List of the centre positions of the surround electrodes.

    """

    # replace electrode name with position if needed
    # and get skin ROI around centre position
    ff = SubjectFiles(fnamehead=fnamehead)
    tmp = sim_struct.ELECTRODE()
    tmp.centre = center_pos
    tmp.substitute_positions_from_cap(ff.get_eeg_cap())

    m = mesh_io.read_msh(fnamehead)
    if tissue_idx < 1000:
        tissue_idx += 1000
    idx = (m.elm.elm_type == 2) & (
        (m.elm.tag1 == tissue_idx) | (m.elm.tag1 == tissue_idx-1000))
    m = m.crop_mesh(elements=m.elm.elm_number[idx])
    P_centre = m.find_closest_element(tmp.centre)
    idx = np.sum((m.nodes[:] - P_centre)**2,
                 1) <= (np.max(radius_surround)+10)**2
    m = m.crop_mesh(nodes=m.nodes.node_number[idx])
    idx = m.elm.connected_components()
    m = m.crop_mesh(elements=max(idx, key=np.size))

    # fit sphere to skin ROI to build local coordinate system
    #   origin: sphere center
    #   x-axis: direction of first surround
    #   z-axis: from sphere center to centre electrode
    r_sph, sph_centre = sphereFit(m.nodes[:])

    M_sph = np.eye(4)
    M_sph[:3, 3] = sph_centre
    tmp = P_centre - sph_centre
    M_sph[:3, 2] = tmp/np.linalg.norm(tmp)
    # direction of first surround
    if pos_dir_1stsurround is not None:
        # replace electrode name with position if needed
        tmp = sim_struct.ELECTRODE()
        tmp.centre = pos_dir_1stsurround
        tmp.substitute_positions_from_cap(ff.get_eeg_cap())
        tmp = tmp.centre - P_centre  # this is not orthogonal to Z
    else:
        # get a vector orthogonal to z-axis
        tmp = np.cross(M_sph[:3, 2], np.eye(3))
        tmp = tmp[:, np.argmax(np.linalg.norm(tmp, axis=1))]
    M_sph[:3, 1] = np.cross(M_sph[:3, 2], tmp)
    M_sph[:3, 1] /= np.linalg.norm(M_sph[:3, 1])
    M_sph[:3, 0] = np.cross(M_sph[:3, 1], M_sph[:3, 2])

    # fit arcs to the skin to get the distances accurate
    if phis_surround is not None:
        if len(phis_surround) != N:
            raise ValueError('exactly N angles are required')
        phis_surround = np.asarray(phis_surround)/180*np.pi  # convert to rad
    else:
        phis_surround = np.arange(N)/N*2*np.pi

    radius_surround = np.array(radius_surround)
    if radius_surround.size == 1:
        radius_surround = np.tile(radius_surround, N)

    N_pts = 50
    P_surround = []
    surround_fit = []
    for phi in range(N):
        theta_on_sph = radius_surround[phi]/r_sph
        arc = np.vstack((r_sph*np.sin(theta_on_sph*np.arange(N_pts)/(N_pts-1)),
                         np.zeros((1, N_pts)),
                         r_sph*np.cos(theta_on_sph *
                                      np.arange(N_pts)/(N_pts-1)),
                         np.ones((1, N_pts))))

        M_rot = np.eye(4)
        M_rot[:3, :3] = R.from_euler('z', phis_surround[phi]).as_dcm()
        M_to_world = M_sph @ M_rot
        M_from_world = np.linalg.inv(M_to_world)

        # project skin points into XZ-plane that contains the arc
        Pts = (M_to_world @ arc).T
        Pts[:, :3] = m.find_closest_element(Pts[:, :3])
        Pts = M_from_world @ Pts.T

        # fit individual arc
        r_arc, arc_centre = sphereFit(Pts[(0, 2), :].T)

        if np.abs(arc_centre[0]) > r_arc:
            # z-axis does not intersect with circle
            # --> use initial sphere instead
            r_arc = r_sph
            arc_centre *= 0

        theta_z0_on_arc = -np.arcsin(arc_centre[0]/r_arc)
        if arc_centre[1] < np.mean(Pts[2, :]):
            theta_on_arc = radius_surround[phi]/r_arc + theta_z0_on_arc
        else:
            # best fitting arc has opposite curvature compared
            # to initial sphere
            theta_z0_on_arc = - theta_z0_on_arc + np.pi
            theta_on_arc = theta_z0_on_arc - radius_surround[phi]/r_arc

        # get centre of surround electrode
        tmp = np.array((r_arc*np.sin(theta_on_arc) + arc_centre[0],
                        0.,
                        r_arc*np.cos(theta_on_arc) + arc_centre[1],
                        1.))
        P_surround.append(m.find_closest_element((M_to_world @ tmp).T[:3]))

        if DEBUG:
            surround_fit.append(
                (arc_centre, r_arc, theta_on_arc, theta_z0_on_arc, M_to_world))
    if DEBUG:
        print('achieved distances:')
        print(np.linalg.norm(np.array(P_surround)-P_centre, axis=1))
        # _show_for_debugging(m, sph_centre, r_sph, P_centre,
        #                     P_surround, surround_fit, M_sph)

    return P_surround


def expand_to_center_surround(S, subpath, radius_surround = 50, N = 4,
                              pos_dir_1stsurround = None, multichannel = False,
                              phis_surround=None, el_surround=None):
    """
    Generate a center-surround montage (standard: 4x1) from a TDCSLIST.

    The TDCSLIST has to contain only the center electrode. Copies of this
    electrode are then placed in a circle around the centre

    Parameters
    ----------
    S : TDCSLIST
        TDCSLIST with the center electrode.
    subpath : string
        m2m_folder of the subject
    radius_surround : float, optional
        Distance (centre-to-centre) between the centre and surround
        electrodes. The default is 50.
    N : int, optional
        Number of surround electrodes. The default is 4.
    pos_dir_1stsurround : array (3,) or string, optional
        A position indicating the direction from center_pos to
        the position of the first surround electrode. The default is None.
    multichannel : Boolean, optional
        When set to True, a multichannel stimulator with each suround channel
        receiving 1/N-th of the of the center channel will be simulated
        (standard: False, i.e. all surround electrodes connected to the
         same return channel).

    Returns
    -------
    S : TDCSLIST
        TDCSLIST with the surround electrodes added.

    """
    if S.type != 'TDCSLIST':
        raise TypeError('The first parameter needs to be a TDCSLIST.')
    if len(S.electrode) != 1:
        raise ValueError('The TDCSLIST has to contain exactly one ELECTRODE.')
    if not os.path.isdir(subpath):
        raise IOError('Could not find m2m-folder: {0}'.format(subpath))

    C = S.electrode[0]
    C.channelnr = 1  # Connect center to channel 1
    if not len(C.name):
        C.name = 'centre'

    # set surround channels and current strengths
    if type(S.currents) == float:
        C_current = S.currents
    else:
        C_current = S.currents[0]

    if multichannel:
        S.currents = -C_current/N*np.ones(N+1)
        S.currents[0] = C_current
        Channel_surround = np.arange(2,N+2)
    else:
        S.currents = [C_current, -C_current]
        Channel_surround = 2*np.ones(N,dtype = int)

    # get centers of surround electrodes
    ff = SubjectFiles(subpath=subpath)
    P_surround = get_surround_pos(C.centre, ff.fnamehead, radius_surround = radius_surround,
                                  N = N, pos_dir_1stsurround = pos_dir_1stsurround,
                                  phis_surround=phis_surround)

    if el_surround is not None:
        C = el_surround

    # get direction vector
    ydir = []
    if len(C.pos_ydir):
        tmp = deepcopy(C)
        tmp.substitute_positions_from_cap(ff.get_eeg_cap())
        ydir = tmp.pos_ydir - tmp.centre

    # add surround electrodes to TDCSLIST
    for i in range(N):
        S.electrode.append(deepcopy(C))
        El = S.electrode[-1]
        El.centre = P_surround[i]
        El.channelnr = Channel_surround[i]
        El.name = 'surround '+str(i+1)
        if len(ydir):
            El.pos_ydir = El.centre + ydir
    return S
# end extra stuff for version 3
####################################


def convert_mask(fn_mask_fsspace, hemi, subpath):
    ''' convert mask roi from fsaverage to individual space and get positions
    '''
    assert hemi in ['lh', 'rh']

    subject_files = SubjectFiles(subpath=subpath)
    if int(__version__[0])>3:
        fn_sphere = get_reference_surf(hemi, 'sphere')
        fn_reg = subject_files.get_surface(hemi, 'sphere_reg')
        fn_central = subject_files.get_surface(hemi, 'central')
    else:
        if hemi == 'lh':
            fn_sphere = templates.cat_lh_sphere_ref
            fn_reg = subject_files.lh_reg
            fn_central = subject_files.lh_midgm
        else:
            fn_sphere = templates.cat_rh_sphere_ref
            fn_reg = subject_files.rh_reg
            fn_central = subject_files.rh_midgm

    surf_sphere = mesh_io.read_gifti_surface(fn_sphere)
    try:
        idx = nib.freesurfer.io.read_label(fn_mask_fsspace)
        idx_mask = np.zeros(surf_sphere.nodes.nr, dtype=np.float32)
        idx_mask[idx] = 1.
    except:
        idx_mask = nib.freesurfer.io.read_morph_data(fn_mask_fsspace)

    idx_mask, _ = transformations._surf2surf(
                idx_mask,
                surf_sphere,
                mesh_io.read_gifti_surface(fn_reg)
                )
    idx_mask = idx_mask > 0.0001
    gm_surf =  mesh_io.read_gifti_surface(fn_central)
    central_pos = gm_surf.nodes[idx_mask]

    return idx_mask, central_pos


def ROI_from_nifti(m, fn_ROI, ROI_shift):
    ''' get a ROI from a nifti file and project on GM surface of mesh
    '''
    gm_label = 1002
    m_GMsurf = m.crop_mesh(tags=gm_label)
    m_GMsurf = m_GMsurf.crop_mesh(elm_type=2)
    m_GMsurf.nodes.node_coord -= ROI_shift*m_GMsurf.nodes_normals().value

    ROI_nifti = nib.load(fn_ROI)
    ROI_image = ROI_nifti.get_fdata() > 0
    ROI_affine = ROI_nifti.get_qform()

    nd = mesh_io.NodeData.from_data_grid(m_GMsurf,ROI_image.astype(np.uint16),
                                         ROI_affine, 'ROI')
    ROI_idx = nd.value > 0
    m_GMsurf.add_node_field(nd,'ROI')
    return m_GMsurf, m_GMsurf.nodes.node_coord[ROI_idx]


def sphericalROI(m, ROI_centerMNI, ROI_radius, ROI_shift, subpath):
    ''' get GM surface inside a spherical ROI
    '''
    gm_label = 1002
    m_GMsurf = m.crop_mesh(tags=gm_label)
    m_GMsurf = m_GMsurf.crop_mesh(elm_type=2)
    m_GMsurf.nodes.node_coord -= ROI_shift*m_GMsurf.nodes_normals().value

    ROI_center = transformations.mni2subject_coords(ROI_centerMNI, subpath)
    # tr_centers = m_GMsurf.elements_baricenters().value
    # ROI_idx = np.linalg.norm(tr_centers-ROI_center,axis=1) <= ROI_radius
    # ed=mesh_io.ElementData(ROI_idx)
    # m_GMsurf.add_element_field(ed,'ROI')
    # return m_GMsurf, tr_centers[ROI_idx]

    ROI_idx = np.linalg.norm(m_GMsurf.nodes.node_coord-ROI_center,axis=1) <= ROI_radius
    nd=mesh_io.NodeData(ROI_idx)
    m_GMsurf.add_node_field(nd,'ROI')
    return m_GMsurf, m_GMsurf.nodes.node_coord[ROI_idx]


def get_closest_skin_pos(pos, m, label_skin = 1005):
    ''' returns the position on the skin that is closest to the
        CoG of the provided positions
    '''
    CoG = np.mean(pos, axis = 0)
    idx_skin = np.where(m.elm.tag1 == label_skin)[0]
    elm_centers = m.elements_baricenters().value[idx_skin]
    distQ=np.sum((elm_centers - CoG)**2,axis=1)
    return elm_centers[np.argmin(distQ)]


# def relabel_internal_air(m, subpath, label_skin = 1005, label_new = 1099,
#                          label_internal_air = 501):
#     ''' relabels skin in internal air cavities to something else;
#         relevant for charm meshes
#     '''
#     subject_files = SubjectFiles(subpath=subpath)
#     # relabel internal skin to some other label
#     label_nifti=nib.load(subject_files.labeling)
#     label_affine=label_nifti.affine
#     label_img = label_nifti.get_fdata().astype(int)
#     label_img = label_img == label_internal_air
#     label_img = mrph.binary_dilation(label_img,iterations=2)

#     m = copy.copy(m)
#     ed = mesh_io.ElementData.from_data_grid(m,label_img,label_affine, order=0)
#     idx = ed.value * (m.elm.tag1 == label_skin)
#     m.elm.tag1[idx] = label_new
#     m.elm.tag2[:] = m.elm.tag1
#     return m


def get_outer_skin_points(m, tol: float = 1e-3, label_skin = 1005):
        """Return indices of points estimated to be on the outer skin surface
        (i.e., not those inside nasal cavities, ear canals etc.). Outer points
        are identified by looking for points which do not intersect the mesh in
        the direction of its normal. This is not perfect but seems to do a
        reasonable job of identifying the relevant points. These may then be
        used for projecting electrodes onto the surface.

        PARAMETERS
        ----------
        tol : float
            Tolerance for avoiding self-intersections.
        label_skin : int
            skin label (standard: 1005)
            
        RETURNS
        -------
        indices : ndarray
            Indices of the outer skin points (0-based).
        """
        assert tol > 0

        skin_faces = m.elm[m.elm.tag1 == label_skin, :3]
        subset = np.unique(skin_faces-1)
        m = mesh_io.Msh(mesh_io.Nodes(m.nodes.node_coord), mesh_io.Elements(skin_faces))

        subset = subset if len(subset) < m.nodes.nr else slice(None)
        n = m.nodes_normals().value[subset]
        # Avoid self-intersections by moving each point slightly along the test
        # direction
        idx = np.unique(m.intersect_ray(m.nodes.node_coord[subset] + tol * n, n)[0][:, 0])
        if isinstance(subset, slice):
            return np.setdiff1d(np.arange(m.nodes.nr), idx, assume_unique=True)
        else:
            return np.setdiff1d(subset, subset[idx], assume_unique=True)


def relabel_internal_air(m, label_skin = 1005, label_new = 1099):
    ''' relabels skin in internal air cavities to something else;
        relevant for charm meshes
    '''
    m = copy.copy(m)
    # outer skin nodes
    idx_skinNodes = get_outer_skin_points(m, label_skin = label_skin) + 1
    
    # internal air triangles
    idx_innerAirTri = (m.elm.elm_type == 2) * (m.elm.tag1 == label_skin)
    idx_innerAirTri *= ~np.any(np.in1d(m.elm.node_number_list, idx_skinNodes).reshape(-1, 4), axis=1)
    
    m.elm.tag1[idx_innerAirTri] = 1099
    m.elm.tag2[:] = m.elm.tag1
    return m


def project_to_pial(pos, m, label_mask = 1100, label_GM = 1002):
    ''' project mask positions to pial surface of tet mesh
    '''
    idx_gm = np.where(m.elm.tag1 == label_GM)[0]
    elm_centers = m.elements_baricenters().value[idx_gm]
    idx_mask_mesh = np.zeros(len(pos), dtype=int)
    for i in range(len(pos)):
        # slow, should be improved at some point
        distQ=np.sum((elm_centers - pos[i])**2,axis=1)
        idx_mask_mesh[i] = np.argmin(distQ)
    idx_mask_mesh = idx_gm[idx_mask_mesh]

    # close small gaps
    idx = np.sum( np.isin(m.elm[idx_gm+1][:,:3],
                          m.elm[idx_mask_mesh+1][:,:3]), axis=1) > 1
    idx_mask_mesh = idx_gm[idx]

    m.elm.tag1[idx_mask_mesh] = label_mask
    m.elm.tag2[:] = m.elm.tag1
    return m


def get_optimal_center_pos(m, label_skin = 1005, label_mask = 1100):
    ''' returns the skin position with the lowest ohmic resistance to the mask
    '''
    # set up and solve FEM
    cond_L = cond.standard_cond()
    cond_list = [c.value for c in cond_L]
    elm_cond = cond.cond2elmdata(m, cond_list)

    S = fem.FEMSystem.tdcs(m, elm_cond, [label_skin, label_mask], [0., 1.],
                           solver_options=None)
    v = S.solve()
    v = mesh_io.NodeData(v, name='v', mesh=m)
    m = fem.calc_fields(v, ['j'], cond=elm_cond)

    # get CoG on skin and project on closest skin position
    cutoff_perc = 0.995 # CoG calculation restricted to skin positions with high J
    idx_skin = np.where(m.elm.tag1 == label_skin)[0]
    elm_sizes = m.elements_volumes_and_areas().value[idx_skin]
    elm_centers = m.elements_baricenters().value[idx_skin]
    magnJ = next(x.value for x in m.elmdata if x.field_name=='magnJ')
    magnJ = magnJ[idx_skin]

    idx_sort = np.argsort(magnJ)
    cum_elm_sizes = np.cumsum(elm_sizes[idx_sort])
    cum_elm_sizes /= cum_elm_sizes[-1]
    idx_cutoff = (cum_elm_sizes >= cutoff_perc).argmax()
    magnJ_cutoff = magnJ[idx_sort][idx_cutoff]

    idx = magnJ >= magnJ_cutoff
    scaleFac = magnJ[idx]*elm_sizes[idx]
    scaleFac = np.reshape(scaleFac,(len(scaleFac),1))
    pos_CoG = np.sum(elm_centers[idx]*scaleFac, axis = 0)/np.sum(scaleFac)

    distQ=np.sum((elm_centers - pos_CoG)**2,axis=1)
    pos_CoG = elm_centers[np.argmin(distQ)]
    return m, pos_CoG


def get_central_gm_with_mask(subpath, hemi, fn_mask_fsspace):
    ''' load bihemispheric GM and add mask as node data
    '''

    subject_files = SubjectFiles(subpath=subpath)

    if int(__version__[0])>3:
        fn_lh_central = subject_files.get_surface('lh', 'central')
        fn_rh_central = subject_files.get_surface('rh', 'central')
    else:
        fn_lh_central = subject_files.lh_midgm
        fn_rh_central = subject_files.rh_midgm

    m_surf = mesh_io.read_gifti_surface(fn_lh_central)
    m_surf.elm.tag1 = 1001 * np.ones(m_surf.elm.nr, dtype=int)
    nr_nodes_lh = m_surf.nodes.nr

    m_rh = mesh_io.read_gifti_surface(fn_rh_central)
    m_rh.elm.tag1 = 1002 * np.ones(m_rh.elm.nr, dtype=int)
    nr_nodes_rh = m_rh.nodes.nr

    m_surf = m_surf.join_mesh(m_rh)
    m_surf.elm.tag2[:] = m_surf.elm.tag1

    idx_mask, _ = convert_mask(fn_mask_fsspace, hemi, subpath)

    if hemi == 'lh':
       idx_mask = np.hstack((idx_mask, np.zeros(nr_nodes_rh,dtype=bool)))
    else:
       idx_mask = np.hstack((np.zeros(nr_nodes_lh,dtype=bool), idx_mask))

    nd=mesh_io.NodeData(idx_mask)
    m_surf.add_node_field(nd,fn_mask_fsspace)
    return m_surf


def run_simus(subpath, pathfem, current_center, N,
              radius_surround, phi_offset,
              EL_center, EL_surround,
              map_to_surf = True, map_to_fsavg = True):
    ''' run a batch of simulations with varying radii and/or phi_offsets
    '''
    S = sim_struct.SESSION()
    S.subpath = subpath
    S.pathfem = pathfem
    S.map_to_surf = map_to_surf
    S.map_to_fsavg = map_to_fsavg
    S.open_in_gmsh = False

    phis_surround = np.arange(N)/N*360
    for radius in radius_surround:
        for offset in phi_offset:
            print('setting up radius '+str(radius)+' phi_offset '+str(offset))
            tdcs_list = S.add_tdcslist()
            tdcs_list.currents = current_center
            center = tdcs_list.add_electrode(copy.copy(EL_center))
            if int(__version__[0])>3:
                tdcs_list.expand_to_center_surround(S.subpath, radius, N,
                                                    pos_dir_1stsurround = center.centre + [0., 0., 20.],
                                                    multichannel=True,
                                                    phis_surround = phis_surround+offset,
                                                    el_surround = EL_surround)
            else:
                tdcs_list = expand_to_center_surround(tdcs_list, S.subpath, radius, N,
                                          pos_dir_1stsurround = center.centre + [0., 0., 20.],
                                          multichannel=True,
                                          phis_surround = phis_surround+offset,
                                          el_surround = EL_surround)
                ff = SubjectFiles(subpath=subpath)
                S.fnamehead = ff.fnamehead
    return S.run()


def analyse_simus(subpath, pathfem, hemi, fn_mask_fsspace,
                  radius_surround, phi_offset, var_name, cutoff):
    ''' get a few key metrics from the results mapped on the central gm
    '''
    res_list = glob(os.path.join(pathfem,'subject_overlays','*.msh'))
    # sort res_list according to simulation number
    numbers = list(list(map(int,re.findall('[0-9]+', x))) for x in res_list)
    idx = np.argsort(np.sum(numbers,axis=1))
    res_list = np.array(res_list)[idx]
    assert len(res_list) == len(radius_surround)*len(phi_offset)

    m_surf = get_central_gm_with_mask(subpath, hemi, fn_mask_fsspace)

    nd_sze = m_surf.nodes_volumes_or_areas().value
    idx_mask = m_surf.nodedata[0].value
    roi_median = np.array([])
    focality = np.array([])
    i=0
    for radius in radius_surround:
        for offset in phi_offset:
            m = mesh_io.read_msh(res_list[i])
            assert m.nodes.nr == m_surf.nodes.nr
            print(os.path.split(res_list[i])[-1])
            nd = next(x.value for x in m.nodedata if x.field_name==var_name)
            m_surf.add_node_field(nd, 'r'+str(radius)+'_p'+str(offset))

            roi_median = np.append(roi_median, np.median(nd[idx_mask]))
            focality = np.append(focality, np.sum(nd_sze[nd > roi_median[-1]]))
            i+=1

    roi_median = roi_median.reshape(len(radius_surround),len(phi_offset))
    focality = focality.reshape(len(radius_surround),len(phi_offset))

    i, j = np.where(roi_median >= cutoff)
    if len(i) == 0:
        raise ValueError("median field strength in ROI never exceeds the threshold!")
    idx = np.argmin(focality[i,j])
    best_condition = [radius_surround[i[idx]], phi_offset[j[idx]]]

    return m_surf, roi_median, focality, best_condition


def analyse_simus2(subpath, pathfem, m_surf,
                   radius_surround, phi_offset, var_name, cutoff):
    ''' get a few key metrics by mapping the results on the given surface
    '''
    res_list = glob(os.path.join(escape(pathfem),'*TDCS*.msh'))
    # sort res_list according to simulation number
    numbers = list(list(map(int,re.findall('[0-9]+', x))) for x in res_list)
    idx = np.argsort(np.sum(numbers,axis=1))
    res_list = np.array(res_list)[idx]
    assert len(res_list) == len(radius_surround)*len(phi_offset)

    var_name, quantity = var_name.split('_',maxsplit = 1)
    if quantity not in ['norm', 'magn', 'normal', 'tangent', 'angle']:
        raise ValueError('Invalid quanty in {0}'.format(quantity))
    print('Analysing '+var_name+' '+quantity)

    def calc_quantities(nd, quantities):
        d = dict.fromkeys(quantities)
        for q in quantities:
            if q == 'norm' or q == 'magn':
                d[q] = nd.norm()
            elif q == 'normal':
                d[q] = nd.normal()
                d[q].value *= -1
            elif q == 'tangent':
                d[q] = nd.tangent()
            elif q == 'angle':
                d[q] = nd.angle()
            else:
                raise ValueError('Invalid quantity: {0}'.format(q))
        return d

    m_surf = copy.deepcopy(m_surf)
    nd_sze = m_surf.nodes_volumes_or_areas().value
    assert len(m_surf.nodedata) == 1
    idx_mask = m_surf.nodedata[0].value > 0
    roi_median = np.array([])
    focality = np.array([])
    i=0
    for radius in radius_surround:
        for offset in phi_offset:
            m = mesh_io.read_msh(res_list[i])
            # Crop out WM, GM, and CSF. We add WM and CSF to make the mesh convex.
            m = m.crop_mesh(tags=[1,2,3])
            # Set the volume to be GM. The interpolation will use only the tetrahedra in the volume.
            th_indices = m.elm.elm_number[m.elm.tag1 == 2]

            # Interpolate to surface
            ed = next(x for x in m.elmdata if x.field_name==var_name)
            assert ed.nr_comp == 3
            interpolated = ed.interpolate_to_surface(m_surf, th_indices=th_indices)
            nd = calc_quantities(interpolated, [quantity]).popitem()[1]
            m_surf.add_node_field(nd, 'r'+str(radius)+'_p'+str(offset))

            roi_median = np.append(roi_median, np.median(nd.value[idx_mask]))
            focality = np.append(focality, np.sum(nd_sze[nd.value > roi_median[-1]]))
            i+=1

    roi_median = roi_median.reshape(len(radius_surround),len(phi_offset))
    focality = focality.reshape(len(radius_surround),len(phi_offset))

    i, j = np.where(roi_median >= cutoff)
    if len(i) == 0:
        raise ValueError("median field strength in ROI never exceeds the threshold!")
    idx = np.argmin(focality[i,j])
    best_condition = [radius_surround[i[idx]], phi_offset[j[idx]]]

    return m_surf, roi_median, focality, best_condition



# def blubb(m_surf):
#     ''' get a few key metrics by mapping the results on the given surface
#     '''
#     m_surf = copy.deepcopy(m_surf)
#     m_surf.add_node_field(m_surf.nodedata[0], 'blubb')
#     m_surf.add_node_field(m_surf.nodedata[0], 'blubb2')
#     return m_surf
