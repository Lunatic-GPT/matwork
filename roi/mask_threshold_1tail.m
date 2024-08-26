function m=mask_threshold_1tail(mag,m,sig_scale,nmax,csize)
% mag: magnitude image
% m: mask for searching the vessel
% sig_scale: the threshold scale
% nmax: maximum number of voxels for the vessel.

  d=mag(m>0);
d=double(d);
  d=sort(d);
if sig_scale>0
  
    mn=mean(d(nmax+1:end));
    sd=std(d(nmax+1:end));
    m=mag>mn+sd*sig_scale & m>0;
else
    mn=mean(d(1:end-nmax));
    sd=std(d(1:end-nmax));
    m=mag<mn+sd*sig_scale & m>0;
end
fprintf('Mean = %f; std = %f\n',mn,sd);
m=clusterize2(m,csize);
m=m>0;


