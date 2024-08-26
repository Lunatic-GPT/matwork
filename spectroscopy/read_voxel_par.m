function [voxc,posc,ntc,time]=read_voxel_par(spc)
%[voxc,posc,ntc,time]=read_voxel_par(spc)
spc=str2cell(spc);
nsc=length(spc);
voxc=zeros(nsc,3);
posc=zeros(nsc,3);
ntc=zeros(1,nsc);
%r=zeros(nsc,11);
%er=zeros(nsc,11);
time = zeros(2,3,nsc);
for i=1:nsc
    voxc(i,1)=readPar(spc{i},'vox1');
    voxc(i,2)=readPar(spc{i},'vox2');
    voxc(i,3)=readPar(spc{i},'vox3');    
    posc(i,1)=readPar(spc{i},'pos1');
    posc(i,2)=readPar(spc{i},'pos2');
    posc(i,3)=readPar(spc{i},'pos3');
    ntc(i)=readPar(spc{i},'nt');
    
    time(:,:,i)=read_time(spc{i});
end
