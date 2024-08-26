function res=swi_calc(im_mag,im_ph,n,sig_smooth)
% swi_calc(im_mag,im_ph,n)
% implement the SWI calculation in E. Mark Haacke et al., Magnetic Resonance in Medicine 52:612–618 (2004)
% Phase unwrap based on Wang et. al., Journal. of Magnetic Resonance Imaging 12:661–670 (2000)

% im_mag: nifti file for magnitude image
% im_ph: nifti file for phase image
% n: number of times to apply the phase mask.  % default n = 1.
% noise_level
if ~exist('n','var')
    n=1;
end


nvox_filter = 64;

if ~isempty(strfind(im_mag,'.nii'))
a=load_untouch_nii(im_ph);
b=load_untouch_nii(im_mag);

out=b;
a=double(a.img);
b=double(b.img);

%%

hdr=load_nii_hdr(im_mag);
voxelsize= hdr.dime.pixdim(2:4);

else
    
    a=ri(im_ph);
    b=ri(im_mag);
end
a=double(a);
b=double(b);
rng=[min(a(:)),max(a(:))];

a=(a-mean(rng))*2*pi/diff(rng);  % convert to [-pi,pi]

im=b.*exp(1i*a);

%im=spatial_smoothing(im,sig_smooth./voxelsize);


if true
fim=fft2c(im);

sz=size(im);
if length(sz)==2
    sz(3)=1;
end

filter=0*fim;

[xx,yy]=meshgrid(0:nvox_filter-1,0:nvox_filter-1);

filter_mat=0.25*(1-cos(2*pi*xx/(nvox_filter-1))).*(1-cos(2*pi*yy/(nvox_filter-1)));


filter(sz(1)/2-nvox_filter/2+1:sz(1)/2+nvox_filter/2,sz(2)/2-nvox_filter/2+1:sz(2)/2+nvox_filter/2,:)=repmat(filter_mat,[1,1,sz(3)]);


im_lowpass=ifft2c(fim.*filter);

ph=angle(im./im_lowpass);

else
    hdr=load_nii_hdr(im_mag);
   voxelsize= hdr.dime.pixdim(2:4);
    [ph Laplacian]=MRPhaseUnwrap(a+pi,'voxelsize',voxelsize);
    ph=ph-pi;
    
end

%%

    mph=(pi-ph)./pi;

  mph(mph>1)=1;  
  res=abs(im).*(mph.^n);

  
if ~isempty(strfind(im_mag,'.nii'))
  out.img=res;
  out.hdr.dime.datatype=4;
  im_mag=strtok(im_mag,'.');
  save_untouch_nii(out,[im_mag,'_swi.nii.gz']);
  
  %%{
  im_ph=strtok(im_ph,'.');
  out.img=ph;
  out.hdr.dime.datatype=16;
  save_untouch_nii(out,[im_ph,'_highpass.nii.gz']);
  
  im_ph=strtok(im_ph,'.');
  out.img=angle(im_lowpass);
  out.hdr.dime.datatype=16;
  save_untouch_nii(out,[im_ph,'_lowpass.nii.gz']);
else
  im_mag=strtok(im_mag,'.');
  save([im_mag,'_swi.mat'],'res');
  
  %%{
  im_ph=strtok(im_ph,'.');
  save([im_ph,'_highpass.mat'],'ph');
  
  im_ph=strtok(im_ph,'.');
  an=angle(im_lowpass);
  save([im_ph,'_lowpass.mat'],'an');
end
    
  
  %}
  %writeanalyze(angle(im_lowpass),[im_ph,'_lowpass']);
  
%writeanalyze(res,[im_mag,'_swi']);    

