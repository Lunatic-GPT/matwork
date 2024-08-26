suma = '/export/res/retino/segmentation/suma_FZong/';
anat = 'anat+orig';
strip_skull = true;

if strip_skull
cmd = sprintf('@SUMA_AlignToExperiment -exp_anat %s -surf_anat %s -strip_skull exp_anat -align_centers',anat,fullfile(suma,'Zong_SurfVol_ns+orig'));
else
    cmd = sprintf('@SUMA_AlignToExperiment -exp_anat %s -surf_anat %s -align_centers',anat,fullfile(suma,'Zong_SurfVol_ns+orig'));
end
%unix(cmd);

p =parameter('convert surface rois to volume');
p = add(p,'string','subject','FZong');
p = add(p,'filename','aligned local anatomy','FZong_SurfVol_ns_Alnd_Exp+orig');

p = add(p,'directoryname','roi directory','rois');
p = add(p,'string','left surface rois','V1.lh.1D.roi');
p = add(p,'string','right surface rois','V1.rh.1D.roi');
%p = add(p,'string','right surface rois','');
p = add(p,'directoryname','suma',suma);

p = add(p,'filename','brain mask','brainmask+orig');

sroi2v_gui_callback(p);
