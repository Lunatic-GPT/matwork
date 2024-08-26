function sp=sp_phasecorr(sp,rp,lp,sw)
% sp_phasecorr(sp,rp[,lp=0])

if ~exist('lp','var')
    lp=0;
    sw=1;
end

w = 0:length(sp)-1;
%theta=rp+lp*w/sw*2*pi;
theta = -rp-lp*w/length(sp);
theta=reshape(theta,size(sp));
sp=sp.*(cos(theta*pi/180)-1i*sin(theta*pi/180));