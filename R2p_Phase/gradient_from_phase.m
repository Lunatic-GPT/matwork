function gr=gradient_from_phase(phase,te,voxelsize)
% gr=gradient_from_phase(phase,te,voxelsize)
% te in units of ms;
% voxelsize in units of mm
% calculate field gradient in each voxel; in units of G/cm
% the output has x,y,z components, corresponding to the gradient in the
% three directions.


voxelsize=voxelsize/10; % convert to cm
te=te/1000; % convert to s;
gamma=4258*2*pi;  % rad/G

sz=size(phase);

gr=zeros([sz(1:3),3]);

tmp=diff(phase,1,1);
gr(2:end-1,:,:,1)=0.5*(tmp(1:end-1,:,:)+tmp(2:end,:,:));
gr([1,end],:,:,1)=tmp([1,end],:,:);

tmp=diff(phase,1,2);
gr(:,2:end-1,:,2)=0.5*(tmp(:,1:end-1,:)+tmp(:,2:end,:));
gr(:,[1,end],:,2)=tmp(:,[1,end],:);

tmp=diff(phase,1,3);
gr(:,:,2:end-1,3)=0.5*(tmp(:,:,1:end-1)+tmp(:,:,2:end));
gr(:,:,[1,end],3)=tmp(:,:,[1,end]);

gr(:,:,:,1)=gr(:,:,:,1)/te/gamma/voxelsize(1);
gr(:,:,:,2)=gr(:,:,:,2)/te/gamma/voxelsize(2);
gr(:,:,:,3)=gr(:,:,:,3)/te/gamma/voxelsize(3);





