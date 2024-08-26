
% PVS
pvsname='18_PVS_rs';
server='andrew';

dcm2nii(pvsname);
winscp_put([pvsname,'.nii.gz']);
t=tic;dcm2pvs(pvsname,server,3);toc(t);

%% acquire the 
nav='';
ScanSliceAdjust(pvsname,nav,'SlicePosition1_base.txt','SlicePosition1.txt');

%%
nav='';
ScanSliceAdjust(pvsname,nav,'SlicePosition2_base.txt','SlicePosition2.txt');


