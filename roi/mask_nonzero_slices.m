function mask_nonzero_slices(fname,dim)
% only works for 3D
ind=1:3;

ind(dim)=[];

a=ri(fname);

b=sum(sum(a,ind(1)),ind(2));

islc=find(b>0);

fprintf('%d slices has non-zero voxels between %d - %d\n',length(islc),islc(1),islc(end));


