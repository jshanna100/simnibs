import pyvista as pv
import numpy as np

def read_curv(fn):
    ''' Reads a freesurfer .curv file
    Parameters
    ------------
    fn: str
        File name
    Returns
    ---------
    curv: np.ndarray
        array with file informatio
    '''
    NEW_VERSION_MAGIC_NUMBER = 16777215

    def read_3byte_integer(f):
        b = f.read(3)
        n = struct.unpack('>i', b'\x00' + b)[0]
        return n

    with open(fn, 'rb') as f:
        # Read magic number as a 3 byte integer
        magic = read_3byte_integer(f)
        if magic == NEW_VERSION_MAGIC_NUMBER:
            vnum = struct.unpack(">i", f.read(4))[0]
            fnum = struct.unpack(">i", f.read(4))[0]
            vals_per_vertex = struct.unpack(">i", f.read(4))[0]
            curv = np.fromfile(f, np.dtype('>f'), vnum)

        else:
            fnum = read_3byte_integer(f)
            curv = np.fromfile(f, np.dtype('>f'), magic) / 100.
    return curv

def elec_plot(mesh, cam_dist=200, return_foc=False):
    # renderer
    plotter = pv.Plotter(off_screen=True)

    # get scalp
    scalp_inds = np.where(mesh.cell_data["gmsh:physical"]==5)[0]
    scalp_cells = mesh.extract_cells(scalp_inds)
    plotter.add_mesh(scalp_cells, color=[.7, .5, .5])

    # get electrodes
    centre_inds = np.where(mesh.cell_data["gmsh:physical"]==101)[0]
    centre_cells = mesh.extract_cells(centre_inds)
    plotter.add_mesh(centre_cells, color="red")

    elec1_inds = np.where(mesh.cell_data["gmsh:physical"]==102)[0]
    elec1_cells = mesh.extract_cells(elec1_inds)
    plotter.add_mesh(elec1_cells, color="blue")

    elec2_inds = np.where(mesh.cell_data["gmsh:physical"]==103)[0]
    elec2_cells = mesh.extract_cells(elec2_inds)
    plotter.add_mesh(elec2_cells, color="blue")

    elec3_inds = np.where(mesh.cell_data["gmsh:physical"]==104)[0]
    elec3_cells = mesh.extract_cells(elec3_inds)
    plotter.add_mesh(elec3_cells, color="blue")

    # camera work
    # centre point of this electrode for focal point
    foc = centre_cells.cell_centers().points.mean(axis=0)
    plotter.camera.focal_point = foc
    # normed vector from origin to focal point
    norm_vec = foc / np.linalg.norm(foc)
    pos = foc + norm_vec * cam_dist
    plotter.camera.position = pos

    image = plotter.screenshot(None, return_img=True)
    if return_foc:
        return image, foc
    else:
        return image

def mag_plot(mesh, foc="elec", cam_dist=200, clim=[0., .6], return_foc=False):
    # renderer
    plotter = pv.Plotter(off_screen=True)

    # get gm
    gm_inds = np.where(mesh.cell_data["gmsh:physical"]==2)[0]
    gm_cells = mesh.extract_cells(gm_inds)
    plotter.add_mesh(gm_cells, scalars="magnE", clim=clim)

    # camera work
    if foc == "elec":
        centre_inds = np.where(mesh.cell_data["gmsh:physical"]==101)[0]
        centre_cells = mesh.extract_cells(centre_inds)
        foc = centre_cells.cell_centers().points.mean(axis=0)

    # centre point of this electrode for focal point
    plotter.camera.focal_point = foc
    # normed vector from origin to focal point
    norm_vec = foc / np.linalg.norm(foc)
    pos = foc + norm_vec * cam_dist
    plotter.camera.position = pos

    image = plotter.screenshot(None, return_img=True)
    if return_foc:
        return image, foc
    else:
        return image

def roi_plot(mesh, foc="centre", cam_dist=200, clim=[0., .6]):
    # renderer
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, show_scalar_bar=False)

    # camera work
    if foc == "centre":
        arr_name = mesh.array_names[0] # assumes first is ROI array
        roi_inds = np.where(mesh[arr_name]==1)[0]
        roi_points = mesh.extract_points(roi_inds)
        foc = roi_points.points.mean(axis=0)
    elif isinstance(foc, (list, tuple, np.ndarray)):
        foc = foc
    else:
        raise ValueError("Unrecognised 'foc' value")

    # centre point of this electrode for focal point
    plotter.camera.focal_point = foc
    # normed vector from origin to focal point
    norm_vec = foc / np.linalg.norm(foc)
    pos = foc + norm_vec * cam_dist
    plotter.camera.position = pos

    image = plotter.screenshot(None, return_img=True)
    return image
