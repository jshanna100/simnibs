 #!/usr/bin/bash

subjs=("001" "003" "004" "005" "006" "008" "009" "012" "013" "016" "017" "018" "020" "022" "024" "025" "027" "028" "029" "030" "032" "033" "034" "035" "062" "064" "066" "072" "074" "079" "080" "081" "086" "089" "091" "095" "097" "098" "105" "106" "994" "995" "996" "997" "998" "999")
subjs=("kocatas" "niemann" "riemann")
orig_dir="/media/Linux5_Data03/hannaj/simnibs/4/"
orig_dir="/home/jev/temp/MeMoSlap/4/"
#orig_dir="/home/jev/simnibs/4/"

for subj in "${subjs[@]}"
do
    echo "$subj"
    cd "$orig_dir""$subj"
    charm $subj T1.nii T2.nii --forcerun --forceqform
done
cd "$orig_dir"
