function mpar=moco_series_motionPar(dname)
% output: n*6: the order of the 6 columns is the same as in the dfile;
% i.e. can be plotted with plot_motion_afni(output);

a=readdPar(dname,'ImageComments',true);

dstr=dir2(dname);
nfile=length(dstr);

mpar=zeros(nfile,6);
for i=2:nfile
    
    b=textscan(a{i},'Motion: %f,%f,%f,%f,%f,%f');
    
    mpar(i,:)=cell2array(b);
    
end

mpar(:,[1,3,5])=-mpar(:,[1,3,5]);

mpar=mpar(:,6:-1:1);

