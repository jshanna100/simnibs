#!/usr/bin/bash

subjs=("001" "003" "004" "005" "006" "008" "009" "012" "013" "016" "017" "018" "020" "022" "024" "025" "027" "028" "029" "030" "032" "033" "034" "035" "062" "064" "066" "072" "074" "079" "080" "081" "086" "089" "091" "095" "097" "098" "105" "106" "994" "995" "996" "997" "998" "999")
#root_dir="/home/jev"
root_dir="/home/hannaj"
orig_dir=$root_dir"/simnibs/3/"
dest_dir=$root_dir"/simnibs/4/"

for subj in "${subjs[@]}"
do
   echo "$subj"
   subj_dest="$dest_dir$subj"
   cp -R "$orig_dir""$subj"/m2m_"$subj"/eeg_positions/ "$subj_dest"/m2m_T1w.nii/
done
