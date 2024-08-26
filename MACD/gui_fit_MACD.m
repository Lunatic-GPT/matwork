function run_fit_MACD(dataName,roiName,dcmDir,par)

if ~exist('par','var')
    par=struct;
end
if ~isfield(par,'FlipAngle')
   par.FA=readdPar(dcmDir,'FlipAngle');
end

if ~isfield(params,'T2')
    par.T2=[28.5,23.5]*0.001;
end

par.roi_vessel=ri(roiName);


p=add(p,'filename','wm mask',''); %leave empty is not needed
p=add(p,'filename','vessel mask','');% save as a mask 
p=add(p,'float','recon vox size',[0.3,0.3]);
p=add(p,'float','acquisition vox size',[0.3,0.3]);
p=add(p,'float','slice thickness',[0.3,0.3]);


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


