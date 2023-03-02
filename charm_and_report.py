# -*- coding: utf-8 -*-
'''
charm section taken over in full from simnibs/cli/charm.py, v.4

'''

import argparse
import os
import shutil
import sys
import textwrap
import time
from simnibs.segmentation import charm_main
from simnibs import __version__, sim_struct, mesh_io
from simnibs.utils.file_finder import SubjectFiles
import Nx1_stuff
import yaml


def parseArguments(argv):

    usage_text = textwrap.dedent('''

CREATE HEAD MESH:
    charm subID T1 {T2}

VISUAL CHECK OF RESULTS:
    open the m2m_{subID}\results.html

RUN ONLY PARTS OF CHARM:
    charm subID T1 T2 --registerT2  (registration of T2 to T1)
    charm subID {T1} --initatlas  (initial affine registration of atlas to MR images)
    charm subID --segment  (make label image, reconstruct surfaces, register to fsaverage and MNI)
    charm subID --mesh  (create head mesh from label images)

    Note: Parts can be concatenated, e.g. charm subID --initatlas --segment

MANUAL EDITING:
    edit m2m_{subID}/label_prep/tissue_labeling_upsampled.nii.gz using a
    viewer of your choice, then call charm subID --mesh to recreate head mesh

    ''')

    parser = argparse.ArgumentParser(
        prog="charm",usage=usage_text)

    parser.add_argument('subID',nargs='?',help="""Subject ID. Charm will create
               the folder m2m_{sub_ID} in the current working directory to
               store the result files.""")
    parser.add_argument('T1',nargs='?',help="T1-weighted image")
    parser.add_argument('T2',nargs='?',help="T2-weighted image (optional)")

    parser.add_argument('-v','--version', action='version', version=__version__)

    parser.add_argument('--registerT2',action='store_true',default=False,
                        help="Register T2- to T1-weighted image")
    parser.add_argument('--initatlas',action='store_true',default=False,
                        help="""Affine registration of atlas to input images
                        (Note:T1-weighted image has to be supplied if no
                        T2-weighted image is used and --registerT2 is thus
                        skipped)""")
    parser.add_argument('--segment',action='store_true',default=False,
                        help="""Run segmentation to create label image,
                        reconstruct the middle cortical surfaces, and create
                        the registrations to the fsaverage and MNI templates""")
    parser.add_argument('--mesh',action='store_true',default=False,
                        help="Create the head mesh from the label image")

    parser.add_argument('--surfaces', action='store_true', default=False,
                        help="Create central cortical surfaces from the label image")

    parser.add_argument('--forcerun',action='store_true',default=False,
                        help="""Overwrite existing m2m_{subID} folder instead
                        of throwing an error""")
    parser.add_argument('--skipregisterT2',action='store_true',default=False,
                        help="""Copy T2-weighted image instead of registering
                        it to the T1-weighted image""")
    parser.add_argument('--usesettings',nargs=1,metavar="settings.ini",
                        help="""ini-file with settings (default: charm.ini in
                        simnibs folder)""")
    parser.add_argument('--noneck', action='store_true', default=False,
                        help="""Inform the segmentation that there's no neck in the scan.""")
    parser.add_argument('--inittransform',  help="""Transformation matrix used
                        to initialize the affine registration of the MNI
                        template to the subject MRI, i.e., it takes the MNI
                        template *to* subject space. Supplied as a path to a
                        space delimited .txt file containing a 4x4
                        transformation matrix (default = None).""")
    parser.add_argument('--forceqform', action='store_true', default=False,
                        help="""Force sform and qform to be the same.""")
    parser.add_argument('--usetransform',  help="""Transformation matrix used
                        instead of doing affine registration of the MNI
                        template to the subject MRI, i.e., it takes the MNI
                        template *to* subject space. Supplied as a path to a
                        space delimited .txt file containing a 4x4
                        transformation matrix (default = None).""")
    parser.add_argument('--debug', action='store_true', default=False,
        help="""Write results from intermediate steps to disk.""")
    args=parser.parse_args(argv)

    # subID is required, otherwise print help and exit (-v and -h handled by parser)
    if args.subID is None:
        parser.print_help()
        exit()

    return args



def main():
    args = parseArguments(sys.argv[1:])
    subject_dir = os.path.join(os.getcwd(), "m2m_"+args.subID)

    if do_charm:
        # run segmentation and meshing

        # check whether it's a fresh run
        fresh_run = args.registerT2
        fresh_run |= args.initatlas and not args.registerT2 and args.T1 is not None # initatlas is the first step in the pipeline when a T1 is explicitly supplied

        if not any([args.registerT2, args.initatlas, args.segment, args.mesh, args.surfaces]):
            # if charm part is not explicitly stated, run all
            fresh_run=True
            args.initatlas=True
            args.segment=True
            args.mesh=True
            args.surfaces = True
            if args.T2 is not None:
                args.registerT2=True

        # T1 name has to be supplied when it's a fresh run
        if fresh_run and args.T1 is None:
            raise RuntimeError("ERROR: Filename of T1-weighted image has to be supplied")

        # T2 name has to be supplied when registerT2==True
        if args.registerT2 and args.T2 is None:
            raise RuntimeError("ERROR: Filename of T2-weighted image has to be supplied")

        if fresh_run and os.path.exists(subject_dir):
            # stop when subject_dir folder exists and it's a fresh run (unless --forcerun is set)
            if not args.forcerun:
                raise RuntimeError("ERROR: --forcerun has to be set to overwrite existing m2m_{subID} folder")
            else:
                if args.usesettings is not None and os.path.dirname(os.path.abspath(args.usesettings[0])) == os.path.abspath(subject_dir):
                    raise RuntimeError("ERROR: move the custom settings file out of the m2m-folder before running with --forcerun.")

                shutil.rmtree(subject_dir)
                time.sleep(2)

        charm_main.run(subject_dir, args.T1, args.T2, args.registerT2,
                       args.initatlas, args.segment, args.surfaces, args.mesh,
                       args.usesettings, args.noneck, args.inittransform,
                       args.usetransform, args.forceqform,
                       " ".join(sys.argv[1:]), args.debug)

    # analyse
    if do_analysis:
        sim_config = yaml.load(args.yaml_path)
        subpath = f"{subject_dir}m2m_{subject_dir}"
        mask, phi, hemi = sim_config["mask"]
        version = int(__version__[0])
        var_name = 'E_magn'

        subject_files = SubjectFiles(subpath=subpath)
        pathfem = os.path.join(subpath, project)
        if os.path.isdir(pathfem):
            raise ValueError(f"{pathfem} already occupied.")
        print(pathfem)
        print(f"Subject {subject_files.subid}")

        m = mesh_io.read_msh(subject_files.fnamehead)
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

        Nx1_stuff.run_simus(subpath, os.path.join(pathfem,'final'),
                            current_center, N, radius, phi,
                            EL_center, EL_surround)

        m_surf, roi_median_f, focality_f, _ = Nx1_stuff.analyse_simus(subpath,
                                                os.path.join(pathfem,'final'),
                                                hemi, mask_path,
                                                [best_radius], [phi],
                                                var_name, cutoff)
        mesh_io.write_msh(m_surf,os.path.join(pathfem,'results_final.msh'))


if __name__ == '__main__':
    main()
