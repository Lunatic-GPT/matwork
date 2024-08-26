function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005


Dx = image([2:end,1],:,:,:) - image;
Dy = image(:,[2:end,1],:,:) - image;
%if size(image,3)>1
%Dz = image(:,:,[2:end,1],:)-image;
%res=cat(5,Dx,Dy,Dz);
%else
res = cat(5,Dx,Dy);
%end


