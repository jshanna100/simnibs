# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 19:56:02 2023

@author: axthi
"""

import os

my_path = os.path.dirname(os.path.realpath(__file__))


class project_template:
    def __init__(self):
        self.proj_nr = 0
        self.roi = ''
        self.hemi = ''              # 'lh', 'rh' or 'cereb'
        self.mask_type = ''         # 'curv' or 'mnivol'
        self.phi = 0.
        self.radius = [0.]

        # --- standard settings shared across projects:
        self.current = 0.002        # current flow through center channel (A)
        self.N_surround = 3         # number of surround electrodes
        self.condition = 'optimal'  # 'closest' or 'optimal'
        # ---

        self.fname_roi = ''         # do not change; auto-generated path to file

    def __setattr__(self, name, value):
        if name == 'radius':
             self.__dict__[name] = [value] if not isinstance(value, list) else value
        else:
             self.__dict__[name] = value

        if name == 'roi':
            self.fname_roi = os.path.join(my_path,'masks',self.roi)

    def __repr__(self):
        return str(vars(self))

    def asdict(self):
        return vars(self)


projects = dict()

# P1
p = project_template()
p.proj_nr = 1
p.roi = 'P1_rTP-RH'
p.hemi = 'rh'
p.mask_type = 'curv'
p.phi = 35.
p.radius = 50.
projects.update({p.proj_nr: p})

# P2
p = project_template()
p.proj_nr = 2
p.roi = 'P2_lPCC-new-LH'
p.hemi = 'lh'
p.mask_type = 'curv'
p.phi = 90.
p.radius = 50.
projects.update({p.proj_nr: p})

# P3
p = project_template()
p.proj_nr = 3
p.roi = 'P3_lTP-LH'
p.hemi = 'lh'
p.mask_type = 'curv'
p.phi = 90.
p.radius = 50.
projects.update({p.proj_nr: p})

# P4
p = project_template()
p.proj_nr = 4
p.roi = 'P4_lIFG-LH'
p.hemi = 'lh'
p.mask_type = 'curv'
p.phi = 75.
p.radius = 50.
projects.update({p.proj_nr: p})

# P5
p = project_template()
p.proj_nr = 5
p.roi = 'P5_lM1-LH'
p.hemi = 'lh'
p.mask_type = 'curv'
p.phi = 90.
p.radius = 50.
projects.update({p.proj_nr: p})

# P6
p = project_template()
p.proj_nr = 6
p.roi = 'P6_cereb.nii.gz'
p.hemi = 'cereb'
p.mask_type = 'mnivol'
p.phi = 90.
p.radius = 50.
projects.update({p.proj_nr: p})

# P7
p = project_template()
p.proj_nr = 7
p.roi = 'P7_rDLPFCnew-RH'
p.hemi = 'rh'
p.mask_type = 'curv'
p.phi = 30.
p.radius = 50.
projects.update({p.proj_nr: p})

# P8
p = project_template()
p.proj_nr = 8
p.roi = 'P8_lDLPFC-LH'
p.hemi = 'lh'
p.mask_type = 'curv'
p.phi = 75.
p.radius = 40.
projects.update({p.proj_nr: p})
