function [residual,cd_fit]=MBAC_fitfunc_mCondition(center, rad,v,fit_par)
% x(1:2) in cm
% x(3) in mm
% x(4) in 10 cm/s

% center: in cm
% rad: in mm
% v in 10 cm/s
% fit_par contains:
    % roi: voxels within roi will be included in fitting
    % data: CD image
    % voxSize: cm
    % voxSize_interp: cm
    % VENC: in cm/s
    % interp_factor: interpolation factor for convolution calculation
    % flow_pattern: Plug, Lamina, or BluntedParobolic
    % heart_rate: min^-1
    % lambda_r
    % lambda_v
    % lambda_rv
lambda_r=fit_par.lambda_r;
lambda_v=fit_par.lambda_v;
lambda_rv=fit_par.lambda_rv;

cd=fit_par.data;
roi= fit_par.roi;
venc=fit_par.VENC;
for i=1:length(rad)
   fit_par.VENC=venc; 
   im1(:,:,1,i)=ComplexImage_PA(center(i,:),rad(i),v(i)*10,fit_par);
   fit_par.VENC=Inf;
   im0(:,:,1,i)=ComplexImage_PA(center(i,:),rad(i),v(i)*10,fit_par);
end
    cd_fit=im1-im0;
    
    resid=getv_roi(cd-cd_fit,roi>0);
%     reg_r=diffc(rad)*lambda_r;
%     reg_v=diffc(v)*lambda_v;
%     reg_rv=lambda_rv*cost_sigmoid(rad,v);   
    residual=[real(resid(:));imag(resid(:))];


%res=[real(resid(:));imag(resid(:))];

function res=diffc(x)

x=x(:);
res=[x(2:end);x(1)]-x;


function res=cost_sigmoid(v,d)

res=sigmoid(-diffc(v).*diffc(d));


function res=sigmoid(x)
res=1/(1+exp(-x));







