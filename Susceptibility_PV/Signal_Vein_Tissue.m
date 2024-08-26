function s=Signal_Vein_Tissue(a,dchi,B0,TE,pos,rhoc,dir_field,dir_vessel)

% B0: unit T
% TE: unit ms
% rhoc is the density at center normalized by rho in tissue
% pos: the position for calculating output; vessel center goes through 0,0;
% the plane normal is defined as z.
% dir_field: 1*3 unit vector for B0
% dir_vessel: 1*3 unit vector for vessel
% 
% a: radius of the vessels; 
% output s: also weighted by rhoc; 


%transform to the cood system where vessel is z and Bo is in the x-z plane
%(phi = 0).

gamma=42.58e6*2*pi;  %rad/T/s
TE = TE/1000;

[theta,phi]=unitVec2thetaPhi(dir_vessel);
xmat=transform_matrix_rotation(theta,phi);

dir_field2=xmat\dir_field(:);
pos2=xmat\pos(:);

[theta_B0,phi]=unitVec2thetaPhi(dir_field2); % theta2 is the angle for B0 wrt vessel
xmat2=transform_matrix_rotation(0,phi);

pos3=xmat2\pos2;
[theta_p,phi_p]=unitVec2thetaPhi(pos3);

r_p = sqrt(pos3(1)^2+pos3(2)^2);

g=0.5*gamma*dchi*B0*TE;

theta_B0=theta_B0*pi/180;

if r_p>a
    phi_p=phi_p*pi/180;
    ps= g*a^2/r_p^2*cos(2*phi_p)*sin(theta_B0)^2;  %make it consistent with Siemens
    s = exp(1i*ps);
    
else
    
    ps=g/3*(3*cos(theta_B0)^2-1); %make it consistent with Siemens
    s = rhoc*exp(1i*ps);
    
end

























