function s=Signal_Vein_Tissue_Fast(a,g,pos,rhoc,xmat,theta_B0,boundary_vein)

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
% boundary_vein:  When r=a, assign to vein (default); otherwise assign to
% tissue



%transform to the cood system where vessel is z and Bo is in the x-z plane
%(phi = 0).
% 
% gamma=42.58e6*2*pi;  %rad/T/s
% TE = TE/1000;
% 
% g=0.5*gamma*dchi*B0*TE;

% [theta,phi]=unitVec2thetaPhi(dir_vessel);
% xmat=transform_matrix_rotation(theta,phi);
% 
% dir_field2=xmat\dir_field(:);
% 
% [theta_B0,phi]=unitVec2thetaPhi(dir_field2); % theta2 is the angle for B0 wrt vessel
% xmat2=transform_matrix_rotation(0,phi);
% 
% pos2=xmat\pos(:);
% 
% pos3=xmat2\pos2;
% 

pos3=xmat\pos;

[theta_p,phi_p]=unitVec2thetaPhi(pos3);

r_p = sqrt(pos3(1,:).^2+pos3(2,:).^2);


theta_B0=theta_B0*pi/180;
phi_p=phi_p*pi/180;
    ps= - g*a^2./r_p.^2.*cos(2*phi_p)*sin(theta_B0)^2;
    
 
    s = exp(1i*ps);
    
     ps=-g/3*(3*cos(theta_B0)^2-1);
     
     ps=-ps; %make it consistent with Siemens
     
    s2 = rhoc*exp(1i*ps);
    
    if boundary_vein
        s(r_p<=a)=s2;
    else
        s(r_p<a)=s2;
    end
    

























