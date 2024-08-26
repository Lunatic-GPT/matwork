function v0=thetaPhi2unitVec(theta,phi,isDeg)

%theta and phi in degree by default;
if ~exist('isDeg','var')
    isDeg=true;
end

if isDeg
 theta=theta*pi/180;
 phi=phi*pi/180;
end

v0=[sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta)];
