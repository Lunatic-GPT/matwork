function p=gui_spv3D()
% perform a three dimensional pattern match

nTE = 2;
p = parameter('SPV3D GUI');
%p = add(p,'filename','mask path file','*.mat'); 
 
p = add(p,'filename','mask file','septalVein.mat'); 
%p = add(p,'button','determine_mask_path','gui_spv_callback(params)');


%p = add(p,'button','calc mom','gui_spv_callback(params)');
% 
% p=add(p,'float','segment length (mm)',0);
% p = add(p,'button','rotate_data','gui_spv_callback(params)');

   p=add(p,'float','B0 (T)',7);
    
%    p=add(p,'filename','Vessel ROI','VesselROI_.mat');
  
%    p=add(p,'filename','Surround ROI','SurroundROI_.mat');

p=add(p,'string','TE for calc susc',num2str(1:nTE));

% 
% p=add(p,'int','TE for calc moment',nTE);

p=add(p,'string','vessel ROI radius (mm)','1.8');

%p=add(p,'float','Mom calc ring radii (mm)',[0.5, 1.2]);

%p=add(p,'float','bg ring radii (mm)',[2,2.5]);

p = add(p,'button','determine_mask_path','gui_spv3D_callback(params)');

p = add(p,'filename','manual mask file','septalVein_rad1.8.mat'); 
p = add(p,'string','vessels',''); % if only fit a certain segment of a vessel; then separate the segment number by a '-'

p = add(p,'button','show pattern','gui_spv3D_callback(params)');
p = add(p,'string','Init. [r,chi]','');
p=add(p,'bool','fit S',false);
p=add(p,'int','interp',5);
p = add(p,'filename','wm mask','roi_wm.mat'); 

p = add(p,'button','calc_from_pattern','gui_spv3D_callback(params)');


% 
 p=add(p,'button','upload data','gui_spv3D_callback(params)');
 p=add(p,'button','upload script','gui_spv3D_callback(params)');
 
 p=add(p,'button','run','gui_spv3D_callback(params)');
 
 p=add(p,'button','upload + run','gui_spv3D_callback(params)');
 
 p=add(p,'button','determine_mask_path+upload+run','gui_spv3D_callback(params)');
 
 p=add(p,'button','get results','gui_spv3D_callback(params)');
 p = add(p,'int','seg',1); % if only fit a certain segment of a vessel; then separate the segment number by a '-'
 p = add(p,'int','1st fig#',100); 

p = add(p,'button','show fit','gui_spv3D_callback(params)');  
%  p=add(p,'filename','saved results','spv_*.mat');
% p = add(p,'button','show results','gui_spv_callback(params)');
p = add(p,'button','save params','gui_spv3D_callback(params)');

p = add(p,'button','Close','close');

if nargout==0
parametergui(p);
p='';
end


