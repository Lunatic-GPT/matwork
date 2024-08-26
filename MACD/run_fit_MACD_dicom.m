function run_fit_MACD_dicom(dcmName,roiName,par)
% par: should at least contain negPhase

cdcmName=str2cell(dcmName);

dcmName=cdcmName{1};

if ~exist('par','var')
    par=struct;
end
if ~isfield(par,'FlipAngle')
   par.FA=readdPar(dcmName,'FlipAngle');
end


if ~isfield(par,'T2')
    par.T2=[28.5,23.5]*0.001;
end

if ~isfield(par,'T1')
    par.T1=[2.6,1.2];
end

if ~isfield(par,'negPhase')
 par.negPhase = false;
end

par.roi_vessel=ri_d1(roiName);

if ~isfield(par,'roi_wm')
   par.roi_wm = ones(size(par.roi_vessel));
end
a=cdp(dcmName);
protocol=[dcmName,'.pro'];


d=get_data(cdcmName);



    par.voxSize=get_acq_voxsize(protocol)*0.1;
    par.voxSize_interp=vec(a.VoxelSize(1:2))'/10; % in cm;
    par.thk=a.Thickness/10;
    par.mag=double(d(:,:,1,2:3));
     
    par.ph = double(d(:,:,1,1))*360/4096;
      
    par.FA=a.FA;
    
    iroiRad=3;
    par.roiRad = par.voxSize_interp(1)*iroiRad;
    par.bgRad = par.voxSize_interp(1)*[11,6];
    par.interp_factor = 32;
    par.VENC = readsPar(protocol, 'sFlowArray.asElm[0].nVelocity');
    
    d_blood=1.06;  %g/ml
    lambda = 1/0.9; %g/ml; blood to brain water partition;
    par.lambda=lambda/d_blood;  %ml/ml;;
    
    par.flow_pattern='BluntedParabolic';

    seg=readsPar(protocol,'lSegments');

    par.TR = a.TR/2/seg/1e3;
    
    par.TE = a.TE/1e3;
   
    par.rad0=0.1;
    par.v0=1;
   
    if ~isfield(par,'savename')
       savename=sprintf('MBAC_%s',dcmName);
    else
       savename=par.savename; 
    end
    Fit_MACD(par,savename);

function data=get_data(fnames)


nscans=length(fnames);

for i=1:nscans
    
    data(:,:,:,:,i)=ri(fnames{i},1);    
end

data=mean(data,5);

% 
% gui_show_flow;
function res=get_acq_voxsize(protocol)
fov(1)=readsPar(protocol,'asSlice[0].dReadoutFOV');
fov(2)=readsPar(protocol,'asSlice[0].dPhaseFOV');

nvox(1)=readsPar(protocol,'lBaseResolution');
nvox(2)=readsPar(protocol,'lPhaseEncodingLines');

res=fov./nvox;