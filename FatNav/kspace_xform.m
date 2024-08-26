
function [k2,d3] = kspace_xform(k,data,xform,voxsize)
% phase encoding is the 2nd dimension
% k: (nro*npe)*3; -0.5 0.5
% data: nro*npe*nch
% xform: 3*4
% fov: 1*3


%k2=xform(:,1:3)*k';
%k2=k2';

k2=k*xform(:,1:3);
k3=reshape(k2,[size(data,1),size(data,2),1,3]);
k3=repmat(k3,[1,1,size(data,3),1]);

shift=xform(:,4)'./voxsize*2*pi; % phase shift
shift=reshape(shift,[1,1,1,3]);

shift=repmat(shift,[size(data),1]);
d3 = data.*exp(-1i*sum(shift.*k3,4));  % 




