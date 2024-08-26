function res=calc_max_motion(an_total,shift,rad)
    
    % total (n*3): 1:2 - theta, phi of axis, 3 - total angle
    % shift (n*3) in mm
    % rad: radius of the sphere in mm;
    
    if ~exist('rad','var')
        
        v=(1094.40+1223.99)/2;  %  https://doi.org/10.1002/hbm.22619
        
        rad=(v/pi*3/4)^(1/3)*10; % 65 mm
    end
    theta=an_total(:,1);
    phi=an_total(:,2);
    
    tot=abs(an_total(:,3));
    
    v0=thetaPhi2unitVec(theta,phi);
    
    z=sum(v0.*shift,2).*v0;
    
    xy=sos(shift-z,2);
    
    xy2=[xy,zeros(length(xy),1)]+rad*[sin(tot*pi/180),cos(tot*pi/180)-1];
    
    xyz=[xy2,z];
    
    res=sos(xyz,2);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    