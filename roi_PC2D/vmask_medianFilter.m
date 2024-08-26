function m=vmask_medianFilter(params)

% params contain: 
% pc: map of PC; assume positive value for vessels; in radians;
% patchSize: size for meidan filter; 
% SNR_mag: magnitude SNR map;
% mask: white matter mask
% alpha_2tail: default 0.05;
% ispc: true - pc data; false - mag data
% thr_factor: 
% use_thr_factor: true or false

if ~isfield(params,'alpha_2tail')
    params.alpha_2tail=0.05;
end


wm=params.wm_mask;
patchSize=params.patchSize;


data=params.data; % in radians;

snr=params.SNR_mag;

if params.ispc
   sd_data=1./snr;
else
   sd_data=double(data)./snr; 
end


data_flt=bg_rm_medianFilter(data,wm,patchSize);

[tmp,sd_data]=bg_rm_medianFilter(sd_data,wm,patchSize); 



z=data_flt./sd_data;

%z=z/std(z(wm>0)); %commented out 5/25/2021
if strcmp(params.method,'threshold_scale')   
 
 m=z>params.thr_scale;
 m(wm==0)=0;
 
elseif strcmp(params.method,'fdr')
 alpha=params.alpha_2tail/2;
 y=cdf('normal',z,0,1);
 p=1-y;
 [~,m_wm]=fdr(p(wm>0),alpha);
 m=wm*0;
 m(wm>0)=m_wm;
 
elseif strcmp(params.method,'Bonferroni')
    interp=params.interp_factor;
    alpha=params.alpha_2tail/2;
    ntest=sum(wm(:)>0)/prod(interp);
    y=icdf('normal',alpha/ntest,0,1);
    m=z>abs(y)*sqrt(2);  %need sqrt(2) because the contrast is a difference between two values, each has STD of sd_data.
    m(wm==0)=0;
end

