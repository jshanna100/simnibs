import nibabel as nib
from nibabel.affines import apply_affine
import matplotlib.pyplot as plt
plt.ion()
import numpy as np

def show_slices(slices):
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap="gray", origin="lower")


suit_coords = np.array([38, -61, -38])
radius = 20
reg_inds = [7, 10]
img = nib.load('/home/jev/Downloads/suit.nii')
vox_coords = apply_affine(np.linalg.inv(img.affine), suit_coords)

img_data = img.get_fdata()

# sphere ROI calculation
# get coordinate array
dims = img_data.shape
y, x, z = np.meshgrid(np.arange(dims[1]),
                      np.arange(dims[0]),
                      np.arange(dims[2]))
xyz = np.stack((x, y, z)).reshape((3, -1))
# distance between coordinates in
dist_vecs = xyz.T - vox_coords
dists = np.linalg.norm(dist_vecs, axis=1)
dists = np.reshape(dists, img_data.shape)
sph_mask = dists < radius

reg_mask = np.zeros_like(img_data, dtype=bool)
for reg_idx in reg_inds:
    this_mask = img_data == reg_idx
    print(this_mask.sum())
    reg_mask += this_mask

final_mask = np.bitwise_and(sph_mask, reg_mask)

mask_img = nib.Nifti1Image(final_mask.astype(np.int8), img.affine)
nib.save(mask_img, "/home/jev/Downloads/mask.nii")

sphere_img = nib.Nifti1Image(sph_mask.astype(np.int8), img.affine)
nib.save(sphere_img, "/home/jev/Downloads/sphere.nii")

reg_img = nib.Nifti1Image(reg_mask.astype(np.int8), img.affine)
nib.save(reg_img, "/home/jev/Downloads/reg.nii")

# slice_0 = anat_img_data[x, :, :]
# slice_1 = anat_img_data[:, y, :]
# slice_2 = anat_img_data[:, :, z]
# show_slices([slice_0, slice_1, slice_2])
