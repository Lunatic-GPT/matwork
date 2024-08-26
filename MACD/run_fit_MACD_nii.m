function run_fit_MACD_nii(fmag,fph,paMask,par)
% fmag should contain 2 images while fph should contain 1.
% the phase unit: 0.01 deg
% par
% par should at least contain .pro

protocol=par.pro;
if ~isfield(par,'FlipAngle')
   par.FA=readsPar(protocol,'adFlipAngleDegree[0]');
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

par.roi_vessel=ri_d1(paMask);
par.roi_vessel=par.roi_vessel(:,:,:,1); %in case there are two identical volumes; to display correctly in ITKSnap

if ~isfield(par,'roi_wm')
   par.roi_wm = ones(size(par.roi_vessel));
end

nii=load_untouch_niigz(fmag);
ph=ri(fph);

par.voxSize=get_acq_voxsize(protocol)*0.1;
par.voxSize_interp=nii.hdr.dime.pixdim(2:3)/10; % in cm;
par.thk= nii.hdr.dime.pixdim(4)/10;
par.mag=double(nii.img);
    
par.ph = double(ph(:,:,1,1))/100;
      
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
    TR_us=readsPar(protocol,'alTR[0]');
    TE_us=readsPar(protocol,'alTE[0]');

    par.TR = TR_us/2/seg/1e6;
    
    par.TE = TE_us/1e6;
   
    par.rad0=0.1;
    par.v0=1;
    nii.img=[];
   par.nii=nii;
   
    if ~isfield(par,'savename')
       savename=sprintf('MBAC_%s',strtok(filename(fmag),'.'));
    else
       savename=par.savename; 
    end
    Fit_MACD(par,savename);

    
% 
% gui_show_flow;
function res=get_acq_voxsize(protocol)
fov(1)=readsPar(protocol,'asSlice[0].dReadoutFOV');
fov(2)=readsPar(protocol,'asSlice[0].dPhaseFOV');

nvox(1)=readsPar(protocol,'lBaseResolution');
nvox(2)=readsPar(protocol,'lPhaseEncodingLines');

res=fov./nvox;