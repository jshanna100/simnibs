#!/bin/bash

component="norm"
subject="001"

outdir="/home/jev/4_emp/"
lh="$outdir/$subject.lh.central"
rh="$outdir/$subject.rh.central"

for simu in "P1", "P2", "P3", "P4", "P5", "P7", "P8"; do


lhres="$outdir/$subject.$simu.$component.lh"
rhres="$outdir/$subject.$simu.$component.rh"


datadir="/home/jev/4_emp/$subject_${simu}/subject_overlays"


###### LEFTVIEW

freeview -f $lh:overlay=$lhres:overlay_threshold=0.05,0.1:overlay_method=linearopaque:overlay_color=colorwheel,inverse:edgethickness=0 $rh:overlay=$rhres:overlay_threshold=0,05,0.1:overlay_method=linearopaque:overlay_color=colorwheel,inverse:edgethickness=0 --screenshot $outdir/${subject}.${simu}.${component} -viewport 3d --camera Azimuth -180 Elevation 180&
done
