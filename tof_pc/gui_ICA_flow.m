function p=gui_ICA_flow()

p = parameter('ICA flow');

p=add(p,'filename','phase file','Phase_MID??.mat');
p=add(p,'filename','mag file','Mag_MID??.mat');
p=add(p,'filename','vessel mask','mask_ICA.mat');

p=add(p,'filename','WM mask','mask_wm.mat');

p=add(p,'filename','protocol','../*.pro');  % for flip angle, vox_size, and VENC
p=add(p,'float','num2deg',0.01);  % for detrend phase data use 1, for reconed data use 0.01;
p=add(p,'float','interp factor',[2,2]);

p=add(p,'float','threshold factor',0.1);  % Minimum intensity (in units of mean intensity in manual ROI) 
% that will be included in the roi

p=add(p,'button','apparent flow','gui_ICA_flow_callback(params)');
p=add(p,'button','-------------undersample recon----------------','');

p=add(p,'float','nonzero fraction','0.2');

p=add(p,'float','lowres interp factor',[2 2]);

p=add(p,'filename','mask for signal removal','');

p=add(p,'filename','recon file','recon_MID??.mat');

p=add(p,'button','make low res','gui_ICA_flow_callback(params)');

p=add(p,'button','group results','edit ICA_flow_highRes_vs_lowRes');



p = add(p,'button','Close','close');

if nargout==0
   parametergui(p);
   p=[];
end  