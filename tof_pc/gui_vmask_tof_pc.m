function p2=gui_vmask_tof_pc()

p = parameter('vessel mask');

p=add(p,'filename','mag file','Mag_MID??_mean.mat');  % leave blank if use the fid file in scan folder.

p=add(p,'filename','phase file','Phase_MID??_mean.mat');

p=add(p,'filename','brain mask','mask_wm_MID??.mat');

p=add(p,'float','brain mask threshold factor',0.02); % mask defined by voxel values greater than 0.05*max


p=add(p,'button','Phase detrend','gui_vmask_tof_pc_callback(params)');
p=add(p,'button','Mag detrend','gui_vmask_tof_pc_callback(params)');
p=add(p,'button','Phase+Mag detrend','gui_vmask_tof_pc_callback(params)');

%% this two files should be located in the folders named as the raw data prefix.
% ok if do not exist.
p=add(p,'bool','Use _detrend',1);
p=add(p,'filename','WM mask','mask_wm_MID??.mat');  % if empty, use brain mask above
p=add(p,'float','phase threshold factor',1.96); % 95% one -tail
p=add(p,'float','mag threshold factor',1.96); % 95% one -tail

p=add(p,'float','max overlap separation (mm)',Inf); % mm2 

p=add(p,'float','shape threshold',5); % 95% one -tail

p=add(p,'float','cluster size thr (mm2)',0.12); % mm2 
p=add(p,'int','max vessel area (mm2)',18); % 95% one -tail

p=add(p,'float','bg circle radius (mm)',1.5); % 95% one -tail

p=add(p,'float','roi circle radius (mm)',0.8); % 95% one -tail

p=add(p,'float','interp factor',[2,3.3333]); % 95% one -tail

p=add(p,'float','voxel size',[0.3125 0.5208]);  % for flip angle, vox_size, and VENC

p=add(p,'filename','protocol','../*.pro');  % for flip angle, vox_size, and VENC

%p=add(p,'float','vesselness threshold (ph)',0.25); % mask defined by voxel values greater than 0.05*max
p=add(p,'float','vesselness threshold (mag)',0.15); % mask defined by voxel values greater than 0.05*max

p=add(p,'int','crop(u d l r)','100 100 200 200');
p=add(p,'button','Frangi','gui_vmask_tof_pc_callback(params)');
p=add(p,'button','Close figures','close(101); close(102); close(103)');

p=add(p,'button','detrend+Frangi','gui_vmask_tof_pc_callback(params)');

p=add(p,'button','int. thr.','gui_vmask_tof_pc_callback(params)');

p=add(p,'button','detrend+int. thr.','gui_vmask_tof_pc_callback(params)');
p = add(p,'button','clear file path','gui_vmask_tof_pc_callback(params)');
p = add(p,'button','Close','close');

if nargout==0
  parametergui(p);
else
    p2=p;
end





