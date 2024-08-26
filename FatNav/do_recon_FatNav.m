function do_recon_FatNav(fname,irep,voxsize,center0)
% irep: the repetitions to process; 1 based;
%% for kernel = 5; rep is 40 is fastest; but memory exceeded; use rep = 10;
%% 9/28/2020: save as nifti file instead of mat file. not finished
% works for Dichter's data; need modifications for other cases;
% center0: whether to set center to 0; should be true for new sequences at VE11; false for Dichter (at least latest ones)
% (PMOCO) except for PVS_R21 and Dichter; default false;
%

if ~exist('center0','var')
    center0=true;
end

while ~exist(fname,'file')
    pause(600);
end


load(fname);


if exist('RepetitionNav','var') %different naming conventions
    Data=navData;
    Line=LineNav;
    Partition=PartitionNav;
    Repetition=RepetitionNav;
    
    
end
pro=strtok(fname,'.');
pro=pro(1:end-7);
pro=name4pat(fullfile(pro,'*.pro'),1);

fov(3)=readsPar(pro,'sCuboid.dThickness');
fov(2)=readsPar(pro,'sCuboid.dPhaseFOV');
fov(1)=readsPar(pro,'sCuboid.dReadoutFOV');
%voxsize=readsPar(pro,'adFree[9]');

norm(1)=readsPar(pro,'sCuboid.sNormal.dSag');
norm(2)=readsPar(pro,'sCuboid.sNormal.dCor');
norm(3)=readsPar(pro,'sCuboid.sNormal.dTra');

pos(1)=readsPar(pro,'sCuboid.sPosition.dSag');
pos(2)=readsPar(pro,'sCuboid.sPosition.dCor');
pos(3)=readsPar(pro,'sCuboid.sPosition.dTra');
dphi=readsPar(pro,'sCuboid.dInPlaneRot');

if center0
    center=[0,0,0];
else
    if norm(1)==1 && dphi==0 %sagital
        center(1)=pos(3)/fov(1);
        center(2)=pos(2)/fov(2);
        center(3)=pos(1)/fov(3);
    elseif norm(1)==1 && abs(dphi-1.57)<0.01
        center(1)=pos(2)/fov(1);
        center(2)=pos(3)/fov(2);
        center(3)=pos(1)/fov(3);
    elseif norm(3)==1 && dphi==0 %transverse
        %center=pos;
       warning('center not set properly');
    elseif norm(3)==1 && abs(dphi-1.57)<0.01
       warning('center not set properly');
        %center=pos([2,1,3]);
    else
        
        error('center cannot be determined');
    end
end

if ~exist('voxsize','var')||isempty(voxsize)
 voxsize=fov(1)/size(Data,1);
end
par.maxLin=round(fov(2)/voxsize)-1; %0 based; without undersampling and partial Fourier
par.maxPar=round(fov(3)/voxsize)-1; %0 based


if ~exist('irep','var') || isempty(irep)
    irep=unique(Repetition)+1;
end



mid=strtok(fname,'.');
%%

par.Lin=Line; %0 based
par.Par=Partition; %0 based
par.Rep=Repetition; %0 based
par.iRep=irep; %1 based
par.iRep_ref=1; %1 based
par.center=center; %taken care in the sequence (at least for VE11)
%(in VB17: No for tse_vfl_pss_FatNav, yes for tse_vfl_pss_FatNav_FNRecon)
par.prefix=[mid,'_recon'];
if ~exist([par.prefix,'.mat'],'file')
  reconFatNav(Data,par);
end

mat2nii_FatNav([par.prefix,'.mat'],pro,voxsize,center0);

%[d,iro] = reconFatNav(navData,LineNav,PartitionNav,RepetitionNav,[mid,'_FatNav'],1:size(navData,1),irep,1,center);


