function p=gui_flfq_ana()

p = parameter('analyze fl fq data');

%p=add(p,'filename','mag file','Mag_MIDxx.mat');  % leave blank if use the fid file in scan folder.

p=add(p,'button','move ph mag','flfq_mvphmag');
p=add(p,'filename','phase file','Phase_MID??.mat');

%% this two files should be located in the folders named as the raw data prefix.
% ok if do not exist.

p=add(p,'filename','WM mask','mask_wm_MID??.mat');


p=add(p,'filename','vessel mask','*_Frangi_vmask.mat');


p=add(p,'filename','protocol','../*.pro');  % for flip angle, vox_size, and VENC

%p =add(p,'float','VENC',4);
 p=add(p,'bool','max pc pixel',1);

%p=add(p,'float','vox size (mm)',[0.3125 0.5]);  % in mm
p=add(p,'float','interp factor',[2,3.3333]);
           
p = add(p,'float','bg size (mm)',3);  % background size estimating baseline phase
p=add(p,'float','num2deg',0.01);  % for detrend phase data use 1, for reconed data use 0.01;
%p = add(p,'int','time point shifts',4);  % shift time points
p = add(p,'int','excluded vessels',0);


p = add(p,'button','clear','gui_flfq_ana_callback(params)');
p=add(p,'int','ind min','1 12');
p=add(p,'int','ind peak','5 6 7');

p=add(p,'int','nshift','5');

p = add(p,'button','Pulsatility time course','gui_flfq_ana_callback(params)');


p=add(p,'string','method','3c');
p=add(p,'float','flip angle scale',1); 

p=add(p,'filename','mean mag file','Mag_MID??_mean.mat');  % leave blank if use the fid file in scan folder.
p=add(p,'filename','mean phase file','Phase_MID??_mean.mat');

%p=add(p,'float','interp factor (4)',1); 
p=add(p,'float','roi radius (mm)',0.9); % crop size before interpolation 
p=add(p,'float','nonzero fraction',0.1);
p=add(p,'int','MID','');
p=add(p,'button','flow','gui_flfq_ana_callback(params)');

p=add(p,'float','intensity adjust',1);

p=add(p,'float','snorm',1);

p=add(p,'button','ICA flow','gui_flfq_ana_callback(params)');
p=add(p,'filename','results','Results_Flow_PartialVolume_method3c.mat');
p=add(p,'button','show fit','gui_flfq_ana_callback(params)');

p=add(p,'filename','mag file','Mag_MID??.mat');

p=add(p,'button','pulsatility (partial volume correction)','gui_flfq_ana_callback(params)');

p=add(p,'button','save params','gui_flfq_ana_callback(params)');

p=add(p,'button','clear file path','gui_flfq_ana_callback(params)');

p = add(p,'button','Close','close');

if nargout==0
   parametergui(p);
   p=[];
   
%gui_vmask_tof_pc;

%gui_resliceT1;
end  

% 
% gui_show_flow;


