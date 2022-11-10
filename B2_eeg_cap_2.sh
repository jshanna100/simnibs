#!/usr/bin/bash

subjs=("001" "003" "004" "005" "006" "008" "009" "012" "013" "016" "017" "018" "020" "022" "024" "025" "027" "028" "029" "030" "032" "033" "034" "035" "062" "064" "066" "072" "074" "079" "080" "081" "086" "089" "091" "095" "097" "098" "105" "106" "994" "995" "996" "997" "998" "999")
#subjs=("001" "005")
dir="/media/Linux5_Data03/hannaj/simnibs"
#dir="/home/jev/simnibs"

for subject in "${subjs[@]}"
do
    simnibs_python electrode_warp_2.py -c $dir/simnibs/EEGcap_incl_cheek_buci.csv -o $dir/4/${subject}/m2m_T1w.nii/eeg_positions/EEGcap_incl_cheek_buci_2 -s $dir/4/${subject}/m2m_T1w.nii
done
cd "$orig_dir"
