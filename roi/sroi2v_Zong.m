function sroi2v_Zong(strip_skull)
% sroi2v_Zong(strip_skull)

suma = '/export/res/retino/suma_zong/';
%suma = '/mnt/windrive/under_shoot/suma_zong/';
if ~exist('strip_skull','var') || strip_skull
strip_skull = true;
anat = 'anat+orig';
else
anat = 'anat_al+orig';    
end


if strip_skull
cmd = sprintf('@SUMA_AlignToExperiment -exp_anat %s -surf_anat %s -strip_skull exp_anat -align_centers',anat,fullfile(suma,'Zong_SurfVol_ns+orig'));
else
    cmd = sprintf('@SUMA_AlignToExperiment -exp_anat %s -surf_anat %s -align_centers -al',anat,fullfile(suma,'Zong_SurfVol_ns+orig'));
end
unix(cmd);

p =parameter('convert surface rois to volume');
p = add(p,'string','subject','Zong');
p = add(p,'filename','aligned local anatomy','Zong_SurfVol_ns_Alnd_Exp+orig');

p = add(p,'directoryname','roi directory',fullfile(suma,'rois022109'));
p = add(p,'string','left surface rois','V1d.lh.1D.roi,V1v.lh.1D.roi,V2d.lh.1D.roi,V2v.lh.1D.roi,V3d.lh.1D.roi,V3vx.lh.1D.roi,V3ax.lh.1D.roi,V4x.lh.1D.roi');
p = add(p,'string','right surface rois','V1d.rh.1D.roi,V1v.rh.1D.roi,V2d.rh.1D.roi,V2v.rh.1D.roi,V3d.rh.1D.roi,V3v.rh.1D.roi,V3ax.rh.1D.roi,V4.rh.1D.roi');
%p = add(p,'string','right surface rois','');
p = add(p,'directoryname','suma',suma);

p = add(p,'filename','brain mask','brainmask+orig');

sroi2v_gui_callback(p);


lh = mask_str('sroi2v_lh+orig');
rh = mask_str('sroi2v_rh+orig');
maskcalc(rh{1:2},lh{1:2},'or(a,b,c,d)',1,'mask_V1');
maskcalc(rh{3:4},lh{3:4},'or(a,b,c,d)',1,'mask_V2');
maskcalc(rh{5:6},lh{5:6},'or(a,b,c,d)',1,'mask_V3');
maskcalc(rh{7},lh{7},'or(a,b)',1,'mask_V3A');
maskcalc(rh{8},lh{8},'or(a,b)',1,'mask_V4');