function gui_spv(nTE)
% perform a three dimensional pattern match

nTE = 2;
p = parameter('SPV GUI');
%p = add(p,'filename','mask path file','*.mat'); 

p = add(p,'filename','mask file','neuriteTracer.mat'); 
%p = add(p,'button','determine_mask_path','gui_spv_callback(params)');


%p = add(p,'button','calc mom','gui_spv_callback(params)');
% 
% p=add(p,'float','segment length (mm)',0);
% p = add(p,'button','rotate_data','gui_spv_callback(params)');

   p=add(p,'float','B0 (T)',7);
    
%    p=add(p,'filename','Vessel ROI','VesselROI_.mat');
  
%    p=add(p,'filename','Surround ROI','SurroundROI_.mat');

p=add(p,'string','TE for calc susc',num2str(1:nTE));

p=add(p,'int','TE for calc moment',nTE);

p=add(p,'string','vessel ROI radius (mm)','1.2');

p=add(p,'float','Mom calc ring radii (mm)',[0.5, 1.2]);

p=add(p,'float','bg ring radii (mm)',[2,2.5]);


p = add(p,'button','calc_from_pattern','gui_spv_callback(params)');

p = add(p,'button','show pattern','gui_spv_callback(params)');
% 
% p=add(p,'button','gen bsub','gui_spv_callback(params)');
%  p=add(p,'filename','saved results','spv_*.mat');
% p = add(p,'button','show results','gui_spv_callback(params)');

p = add(p,'button','Close','close');


parametergui(p);


