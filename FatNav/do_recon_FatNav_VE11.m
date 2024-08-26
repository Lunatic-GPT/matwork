function do_recon_FatNav_VE11(fname,maxLin,maxPar,irep)
% irep: the repetitions to process; 1 based;
%% for kernel = 5; rep is 40 is fastest; but memory exceeded; use rep = 10;
%% 9/28/2020: save as nifti file instead of mat file. not finished
%fname: *_FatNav.mat file from SiemensRawData2mat.m
%maxLin: image size -1 along PE
%maxPar: image size - 1 along PAR
%irep: the repetitions to recon; 1 based; default all;


load(fname);

%1 based in VE11, convert to 0 based
Repetition=Repetition-1;
Line=Line-1;
Partition=Partition-1;

if ~exist('irep','var')
    irep=unique(Repetition)+1;
end



mid=strtok_no(fname,'_',2);    
%%

par.Lin=Line; %0 based
par.Par=Partition; %0 based
par.Rep=Repetition; %0 based
par.iRep=irep; %1 based
par.iRep_ref=1; %1 based
par.center=[0,0,0]; %taken care in the sequence (at least for VE11)
par.prefix=[mid,'_FatNav'];
par.maxLin=maxLin; %0 based; without undersampling and partial Fourier 
par.maxPar=maxPar; %0 based

reconFatNav(Data,par);  

%[d,iro] = reconFatNav(navData,LineNav,PartitionNav,RepetitionNav,[mid,'_FatNav'],1:size(navData,1),irep,1,center);  


