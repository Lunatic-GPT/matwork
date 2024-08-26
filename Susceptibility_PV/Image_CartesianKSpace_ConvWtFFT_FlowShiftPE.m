function [res,res0]=Image_CartesianKSpace_ConvWtFFT_FlowShiftPE(fr,fr2,FOV,voxSize,voxSize_interp,interp)
% res=Image_CartesianKSpace_ConvWtFFT(fr,center,FOV,voxSize,voxSize_interp)
% fr: the function of f(x,y) or f(x,y,z); 
% x,y,z are matrices with the same size. x, y, z are the first, second, and
% third dimensions, respectively.
% FOV: field of view; mm; the center of FOV is assumed to be [0,0,0];
% voxSize: acquired voxel size; [2] or [3]
% voxSize_interp: voxel size for output res;
% interp: interpolation factor for calculating the PSF; relative to
% voxSize_interp
% fr2: flow induced shift 
if ~exist('interp','var')
    interp=4;
end

nd=length(voxSize);

if length(FOV)==1
    FOV=FOV*ones(1,nd);
end


if length(voxSize_interp)==1
    voxSize_interp=voxSize_interp*ones(1,nd);
end

%FOV = single(FOV);
nacq=FOV./voxSize;

nvox=round(nacq.*voxSize./voxSize_interp);

%kmax=dk.*nacq/2;
kmax=pi./voxSize;
kmax_interp=kmax.*voxSize./voxSize_interp.*interp;


%dxy_interp=dxy/interp;

% xk=linspace(-lim(1),lim(1),round(lim(1)/dxy(1)*5));
% yk=linspace(-lim(2),lim(2),round(lim(2)/dxy(2)*5));

x0=linspace(-FOV(1)/2,FOV(1)/2,nvox(1)*interp+1);
y0=linspace(-FOV(2)/2,FOV(2)/2,nvox(2)*interp+1);

x0=x0(1:end-1);
y0=y0(1:end-1);

if nd==3
    z0=linspace(-FOV(3)/2,FOV(3)/2,nvox(3)*interp+1);
    z0=z0(1:end-1);
    nz=length(z0);
    z0=reshape(z0,[1,1,length(z0)]);
    z=repmat(z0,[length(x0),length(y0),1]);
else
    nz=1;
end

x=repmat(x0',[1,length(y0),nz]);
y=repmat(y0,[length(x0),1,nz]);
% 

if nd==3
   res0=fr(x,y,z);
   fres0=fft2c(res0);
   fres0=fft1c(fres0,3);
else
   res0=fr(x,y); 
  % fres0=fft2c(res0);
   shift=fr2(x,y);
           
   sz=size(x);
   l=(1:size(x,2))-size(x,2)/2-1;
   fres0=0*res0;
   for i=1:size(x,1)
       for j=1:size(x,2)
           
           if res0(i,j)==0
               continue;
           end
         
           fres0tmp=res0(i,j)*exp(-1i*2*pi*l'*(i-sz(1)/2-1)/sz(1))*exp(-1i*2*pi*l*(j-sz(2)/2-1)/sz(2))/sqrt(sz(1)*sz(2));
         
           
%            tic;
%            res0tmp=res0*0;
%            res0tmp(i,j)=res0(i,j);
%            fres0tmp=fft2c(res0tmp);
%            toc;

          if datenum(version('-date'))>=736580
           fres0=fres0+fres0tmp.*exp(-1i*shift(i,j)/FOV(2)*2*pi*l);
          else
              fres0=fres0+fres0tmp.*exp(-1i*shift(i,j)/FOV(2)*2*pi*repmat(l,[length(l),1]));
          end
       end
   end
   
end

kx=linspace(-kmax_interp(1),kmax_interp(1),nvox(1)*interp+1);
ky=linspace(-kmax_interp(2),kmax_interp(2),nvox(2)*interp+1);
kx=kx(1:end-1);
ky=ky(1:end-1);
if nd==3
    kz=linspace(-kmax_interp(3),kmax_interp(3),nvox(3)*interp+1);
    kz=kz(1:end-1);
end

if nd==3
    fres0(kx<-kmax(1) | kx>=kmax(1),:,:)=0;
    fres0(:,ky<-kmax(2) | ky>=kmax(2),:)=0;
    fres0(:,:,kz<-kmax(3) | kz>=kmax(3))=0;
    
    res=ifft2c(fres0);
    res=ifft1c(res,3);
    res=res(1:interp:end,1:interp:end,1:interp:end);
else
    fres0(kx<-kmax(1) | kx>=kmax(1),:)=0;
    fres0(:,ky<-kmax(2) | ky>=kmax(2))=0;
    res=ifft2c(fres0);
    res=res(1:interp:end,1:interp:end);
end




