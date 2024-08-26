function kdata=kspaceData_circleImage(fr,angular_dependent,nvox,radius,kx,ky)
% nvox: number of voxel within the circle
% res=circleImage_CartesianKSpace(fr,center,FOV,voxSize,voxSize_interp)
% fr: the function of f(radius) or f(radius,phi) where phi is measured from
% x axis which is the second dimension.
% center: the center of the circle; mm
% FOV: field of view; mm
% voxSize: acquired voxel size
% voxSize_interp: voxel size for output res;
% interp: interpolation factor for calculating the PSF; relative to
% voxSize_interp
% 

if ~exist('angular_dependent','var')
    angular_dependent = false;
end

if ~exist('interp','var')
    interp=10;
end

x0=linspace(-radius,radius,nvox);
y0=x0;
x=repmat(x0,[length(y0),1]);
y=repmat(y0',[1,length(x0)]);

res0=fr2(fr,x,y,angular_dependent);

kdata=zeros(size(kx,1),size(kx,2));
for i=1:size(kx,1)
    disp(i);
    for j=1:size(kx,2)  
        tmp=res0.*exp(-1i*kx(i,j)*x).*exp(-1i*ky(i,j)*y);
        kdata(i,j)=mean(tmp(:));
    end
end

function res=fr2(fr,x,y,angular_dependent)

if angular_dependent
   phi= atan2(y,x)*180/pi;
   r=sqrt(x.^2+y.^2);
   res=fr(r,phi);
    
else
    
r=sqrt(x.^2+y.^2);
res=fr(r);
end
% 
% function res=fr2(rad,x,y)
% 
% res=ones(size(x));
% r=sqrt(x.^2+y.^2);
% res(r>rad)=0;




