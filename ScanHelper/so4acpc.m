function so4acpc(acpc,shift)
% shift is the shift (in mm) in the superior direction relatiave to the
% acpc scan
 
if ~exist('shift','var')
    shift=15;
end

ori=readdPar(acpc,'ImageOrientationPatient');
sroi.vec_lr=ori(1:3);
norm=cross(ori(1:3),ori(4:6));
[~,sroi.center]=dcmDimCenter(acpc);

if norm(3)<0
    norm=-norm;
end
sroi.norm=norm;


save_slicePosition(1,sroi.center+shift*sroi.norm,sroi.norm);


function save_slicePosition(wip,pos,norm)

if exist('Y:\','dir')
 fname=sprintf('Y:\SlicePosition%d.txt',wip);
else
 fname=sprintf('SlicePosition%d.txt',wip);
end

fid=fopen(fname,'w');

fprintf(fid,'position=(%f,%f,%f)\n',pos);
fprintf(fid,'normal=(%f,%f,%f)\n',norm);
fclose(fid);

copyfile(fname,filename_append(fname,'_ACPC'));
