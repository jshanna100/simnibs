import os
import numpy as np
import nibabel as nib
from simnibs import SIMNIBSDIR
from simnibs.utils.transformations import volumetric_nonlinear
from simnibs.utils.transformations import volumetric_affine
from simnibs.utils.settings_reader import read_ini
from simnibs.mesh_tools.meshing import create_mesh
from simnibs.mesh_tools.mesh_io import write_msh


def _map_rois(roi_path, sub_path, rois):

    target_scan = os.path.join(sub_path, "T1.nii.gz")
    target_im = nib.load(target_scan)
    target_dim = list(target_im.get_fdata().shape)
    im_deformation = nib.load(
        os.path.join(sub_path, "toMNI", "Conform2MNI_nonl.nii.gz")
    )
    def_tuple = (im_deformation.get_fdata(), im_deformation.affine)

    for i, b in enumerate(rois):
        print("Processing atlas: " + b)
        im_atlas = nib.load(os.path.join(roi_path, b))
        transformed_data = volumetric_nonlinear(
            (im_atlas.get_fdata(), im_atlas.affine), def_tuple, intorder=0
        )
        tranformed_data = np.squeeze(transformed_data)

    return transformed_data


seg_path = "/home/jev/simnibs/ROI/"
subject_path = "/home/jev/simnibs/4/005/m2m_005/"

rois = [
    "P6.nii.gz",
]

roi_seg = _map_rois(seg_path, subject_path, rois)
# Write the segmentation to disk
im_tmp = nib.load(os.path.join(subject_path, "T1.nii.gz"))
im_seg = nib.Nifti1Image(roi_seg.astype(float), im_tmp.affine)
nib.save(im_seg, os.path.join(subject_path, "roi_segmentation.nii.gz"))

# Add the bundle seg to the upsampled tissue seg
im_upsampled_tissue_seg = nib.load(
    os.path.join(subject_path, "label_prep", "tissue_labeling_upsampled.nii.gz")
)
tissue_seg = im_upsampled_tissue_seg.get_fdata()
upsampled_bundles = volumetric_affine(
    (roi_seg, im_tmp.affine),
    np.eye(4),
    im_upsampled_tissue_seg.affine,
    tissue_seg.shape,
    intorder=0,
)

roi_labels = np.unique(roi_seg)
roi_labels = roi_labels[roi_labels > 0]

for label in roi_labels.tolist():
    tissue_seg[upsampled_bundles == label] = 50 + label

im_upsampled_added = nib.Nifti1Image(tissue_seg, im_upsampled_tissue_seg.affine)
nib.save(
    im_upsampled_added,
    os.path.join(
        subject_path, "label_prep", "tissue_labeling_upsampled_with_roi.nii.gz"
    ),
)

# Remesh the new tissue_labels
settings = read_ini(os.path.join(SIMNIBSDIR, "charm.ini"))["mesh"]

new_mesh = create_mesh(
    tissue_seg.astype(np.uint16),
    im_upsampled_tissue_seg.affine,
    elem_sizes=settings["elem_sizes"],
    smooth_size_field=settings["smooth_size_field"],
    skin_facet_size=settings["skin_facet_size"],
    facet_distances=settings["facet_distances"],
    optimize=settings["optimize"],
    remove_spikes=settings["remove_spikes"],
    skin_tag=settings["skin_tag"],
    hierarchy=settings["hierarchy"],
    smooth_steps=settings["smooth_steps"],
    sizing_field=None,
)

write_msh(new_mesh, os.path.join(subject_path, "mesh_with_cereb_roi.msh"))
