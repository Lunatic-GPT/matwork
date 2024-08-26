function [DCF,gdata2]=calc_DCF(crds,effMtx)
% crds: a 3*N*n matrix specifying the cordinate of the sampled k-space locations. 
% reference: James Pipe et al., Magnetic Resonance in Medicine 41:179–186 (1999)

sz=size(crds);
gdata=grid3_MAT_samegrid_real(ones(sz(2:end)),crds,effMtx);
gdata=1./gdata;
for i=1:3

gdata2=grid3_MAT_samegrid_real(gdata,crds,effMtx);

gdata=gdata./(gdata2).^2;
end

DCF=gdata;