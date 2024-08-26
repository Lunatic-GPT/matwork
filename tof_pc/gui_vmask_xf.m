function gui_vmask_xf()

p = parameter('gui to generate vessel mask');

p=add(p,'filename','mag file','*nocad_80');  % leave blank if use the fid file in scan folder.

p=add(p,'filename','phase file','*nocad_80_ph');

p=add(p,'filename','brain mask','mask_wm.mat');
p=add(p,'bool','mask: flip xy',1);  

%p=add(p,'float','brain mask threshold factor',0.02); % mask defined by voxel values greater than 0.05*max
p=add(p,'button','Phase detrend','gui_vmask_tof_pc_callback(params)');
p=add(p,'button','Mag detrend','gui_vmask_tof_pc_callback(params)');
p=add(p,'button','Phase+Mag detrend','gui_vmask_tof_pc_callback(params)');

%% this two files should be located in the folders named as the raw data prefix.
% ok if do not exist.
p=add(p,'bool','Use _detrend',1);
%p=add(p,'filename','WM mask','');  % if empty, use brain mask above
p=add(p,'float','phase threshold factor',-1.96); % 95% one -tail
p=add(p,'float','mag threshold factor',1.96); % 95% one -tail

p=add(p,'float','max overlap separation (mm)',0.6); % mm2 

p=add(p,'float','shape threshold',2); % 95% one -tail

p=add(p,'float','cluster size thr (mm2)',0.12); % mm2 
p=add(p,'int','max vessel area (mm2)',18); % 95% one -tail


p=add(p,'float','bg circle radius (mm)',1.5); % 95% one -tail

p=add(p,'float','roi circle radius (mm)',0.8); % 95% one -tail

%p=add(p,'float','interp factor',[1,1]); % 95% one -tail

p=add(p,'float','voxel size',[0.1953 0.1953]);  % for flip angle, vox_size, and VENC

%p=add(p,'filename','protocol','');  % for flip angle, vox_size, and VENC

p=add(p,'button','gen vessel mask','gui_vmask_tof_pc_callback(params)');

p=add(p,'button','detrend+gen vessel mask','gui_vmask_tof_pc_callback(params)');

p = add(p,'button','clear file path','gui_vmask_tof_pc_callback(params)');

p=add(p,'filename','vessel mask','');
p=add(p,'float','VENC',-4);
p=add(p,'button','calc flow','gui_vmask_tof_pc_callback(params)');

p = add(p,'button','Close','close');

p = parametergui(p);






