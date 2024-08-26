function r=sp_harm(theta,phi)
%r=sp_harm(name,theta,phi)
% the coefficients of components
% {'x1','y1','z1','xz','yz','z2','xy','x2y2'}; in the imaging cordinate
% system, which is different for the shim coil system.
r=zeros(1,8);
    r(1)=sin(theta)*cos(phi);
    r(2)=sin(theta)*sin(phi);
    r(3)=cos(theta);
    r(4)=sin(theta)*cos(theta)*cos(phi);
    r(5)=sin(theta)*cos(theta)*sin(phi);
    r(6) = (3*cos(theta)*cos(theta)-1)/2;
    r(7)=sin(theta)*sin(theta)*sin(2*phi);
    r(8)=sin(theta)*sin(theta)*cos(2*phi);
    
    
    
    
    
