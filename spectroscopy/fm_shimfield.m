function field=fm_shimfield(theta,phi,k,pos,dr)
% field=fm_shimfield(theta,phi,k,pos,dr)
% k are the coefficients of the components {'x1','y1','z1','xz','yz','z2','xy','x2y2'};

r=pos+dr*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
r=sqrt(r(1).^2+r(2).^2+r(3).^2);
mr=[r,r,r,r*r,r*r,r*r,r*r,r*r];

k=reshape(k,[1,8]);
field = sum(k.*sp_harm(theta,phi).*mr);






