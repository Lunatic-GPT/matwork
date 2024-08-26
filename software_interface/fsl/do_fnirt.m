bet $1 $1_bet
flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in $1_bet

function do_fnirt(fin)


cmd=sprintf('!fnirt --ref=MNI152_T1_2mm.nii --in=my_brain.nii --subsamp=4,2,1 ----reffwhm=8,4,0',