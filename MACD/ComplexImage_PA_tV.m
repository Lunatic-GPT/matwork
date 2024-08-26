function im=ComplexImage_PA_tV(center,rad,Vt,fit_par)
% center in cm; 1*2 or n*2
% rad in mm
% vmean in cm/s
% fit_par contains:
% roi: voxels within roi will be included in fitting
% sz: matrix size of the image
% voxSize: cm
% voxSize_interp: cm
% VENC: in cm/s
% interp_factor: interpolation factor for convolution calculation
% va: array of blood velocities cm/s
% sa: signal intensity at va 

if size(center,2)==1
    center=center(:)';
end

if size(center,1)==1
    center=repmat(center,[length(Vt),1]);
end

rad=rad*0.1; % to cm

sz = fit_par.sz;
voxSize=fit_par.voxSize;
voxSize_interp=fit_par.voxSize_interp;
venc=fit_par.VENC;
flow_pattern=fit_par.flow_pattern;

va=fit_par.va;  %cm/s
sa=fit_par.sa;
dvdt=fit_par.dvdt;


FOV=sz.*voxSize_interp;

Tc=60/fit_par.heart_rate;


for i=1:length(Vt)
    dvdt_tmp=get_dvdt(Vt,Tc,i);    
    fr2=@(r) fr(r,Vt(i),dvdt_tmp,va,dvdt,sa,rad(i),venc,flow_pattern);
    im(:,:,:,i)=circleImage_CartesianKSpace(fr2,center(i,:),FOV,voxSize,voxSize_interp,false,fit_par.interp_factor);
end

function dvdt=get_dvdt(Vt,Tc,i)

  n=length(Vt);
  i2=mod(i-2,n)+1;
  dvdt=(Vt(i)-Vt(i2))/Tc*n;
%   if dvdt>2
%       dvdt=2;
%   elseif dvdt<-2
%       dvdt=-2;
%   end
%  
  
function res=fr(r,v,dvdt_tmp,va,dvdt,sa,rad,venc,flow_pattern)
if strcmp(flow_pattern,'Laminar')
    vi=v_laminarFlow(r/rad,v);
elseif strcmp(flow_pattern,'BluntedParabolic')
    vi=v_BluntedParabolicFlow(r/rad,v);
elseif  strcmp(flow_pattern,'Plug')
    vi=v*ones(size(r));
else
    error('Unknow velocity pattern');
end

[va,dvdt]=meshgrid(va,dvdt);
%warning('check here');
dvdt_i=vi/v*dvdt_tmp;
dvdt_i(dvdt_i>2)=2;
dvdt_i(dvdt_i<-2)=-2;

si=interp2(va,dvdt,sa,vi,dvdt_i);
tmp=interp2(va,dvdt,sa,vi,dvdt_i*0);
si(isnan(si))=tmp(isnan(si));
res=si.*(exp(1i*vi/venc*pi));
res(r>rad)=0;

if any(isnan(res(:)) | isinf(res(:)))
    fprintf('dvdt_tmp = %f; rad = %f\n',dvdt_tmp,rad);
end
    
    

