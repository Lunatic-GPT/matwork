function res=swi_calc_Kelly(im_mag,im_ph,n,save_inv)
% swi_calc(im_mag,im_ph,n)
% implement the SWI calculation in E. Mark Haacke et al., Magnetic Resonance in Medicine 52:612–618 (2004)
% Phase unwrap based on Wang et. al., Journal. of Magnetic Resonance Imaging 12:661–670 (2000)

% im_mag: nifti file for magnitude image
% im_ph: nifti file for phase image
% n: number of times to apply the phase mask.  % default n = 1.
% whether to save inverted images


if ~exist('n','var')
    n=1;
end 

if ~exist('save_inv','var')
    save_inv=false;
end


nvox_filter = 64;

a=load_nii(im_ph);
b=load_nii(im_mag);

out=b;
a=double(a.img);
b=double(b.img);

%%


rng=[min(a(:)),max(a(:))];

a=(a-mean(rng))*2*pi/diff(rng);  % convert to [-pi,pi]

im=b.*exp(1i*a);

%im=spatial_smoothing(im,sig_smooth./voxelsize);

%% Remove background phase 
fim=fft2c(im);

sz=size(im);

filter=0*fim;

[xx,yy]=meshgrid(0:nvox_filter-1,0:nvox_filter-1);

filter_mat=0.25*(1-cos(2*pi*xx/(nvox_filter-1))).*(1-cos(2*pi*yy/(nvox_filter-1)));


filter(sz(1)/2-nvox_filter/2+1:sz(1)/2+nvox_filter/2,sz(2)/2-nvox_filter/2+1:sz(2)/2+nvox_filter/2,:)=repmat(filter_mat,[1,1,sz(3)]);


im_lowpass=ifft2c(fim.*filter);

ph=angle(im./im_lowpass);

%% Apply phase mask n times.

  mph=(pi-ph)./pi;
  mph(mph>1)=1;  
  res=b.*(mph.^n);

  
  out.img=res;
  out.hdr.dime.datatype=4;
  im_mag=strtok(im_mag,'.');
  save_nii(out,[im_mag,'_swi.nii']);
  
  if save_inv
     
      out.img=max(res(:))-res;
      save_nii(out,[im_mag,'_swi_inv.nii']);
      
  end
  
  %}
  %writeanalyze(angle(im_lowpass),[im_ph,'_lowpass']);
  
%writeanalyze(res,[im_mag,'_swi']);    

function res = ifft2c(x)

S = size(x);
fctr = S(1)*S(2);

x = reshape(x,S(1),S(2),prod(S(3:end)));

res = zeros(size(x));
for n=1:size(x,3)
res(:,:,n) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n))));
end


res = reshape(res,S);

function res = fft2c(x)

S = size(x);
fctr = S(1)*S(2);

x = reshape(x,S(1),S(2),prod(S(3:end)));

res = zeros(size(x));
for n=1:size(x,3)
	res(:,:,n) = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x(:,:,n))));
end

res = reshape(res,S);
