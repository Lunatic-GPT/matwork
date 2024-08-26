function [res,flag,resid_magMoment]= calc_vessel_susceptibility_nullVein(cd,mv,R,theta_B0,TE,B0,vox_size)
% Note when vein signal is small, the result is not reliable; strongly
% depend on initial condition.
% assume B0 projection is along the first dimension 
% angles in theta_B0 is in degree
% angle in d is in radian
% TE in ms;
% theta_B0 in degree
% vox_size in mm, one element; assume isotropic;
% res: 1*3; res(1) is dchi (ppm); res(2) is radius in mm; res(3) is dchi*radius^2 (unit: ppm*mm2).
% mv: vessel mask
% cd: complex data
% R: radius in mm
% 
% 
% mg=ri(d,[],[],'mg');
% p=ri(d,[],[],'p');
% 
% p=permute(p,[2,1,3]);
% mg=permute(mg,[2,1,3]);
% 
% cd=mg.*exp(1i*p);

cd=permute(cd,[2,1,3]);
p=angle(cd);

mv=permute(mv,[2,1,3]);

y=zeros(1,length(R));

m_ind=ind2subb(size(mv),find(mv>0));
center=mean(m_ind,1);

for i=1:length(R)

    m=mask_circle(size(cd),R(i)/vox_size,center,1);
    y(i)=mean_roi(cd,m)*sum(m(:)>0);
           
end
 

ml=mask_circle(size(cd),R(1)/vox_size,center,1);
mo=ml>0&mv==0;

[dchi_a2,resid_magMoment]=magMoment_outSidePattern(p*180/pi,TE,B0,theta_B0,mo>0,center,vox_size*ones(1,3));

%func=@(b) TotalSignal_Tissue_PerpPlane(b(2),dchi_a2/b(2)/b(2),B0,TE,theta_B0,R)*b(1)-real(y);
%[tmp,t2,t3,flag]=lsqnonlin(func,[rho0,0.4],[0,0],[6*rho0,0.4]);
%    res(1)=dchi_a2/tmp(2)/tmp(2);
%    res(2) = tmp(2);  %in units of mm 

r0=0.2;
rho0=real(y(1))/TotalSignal_Tissue_PerpPlane(r0,dchi_a2/r0/r0,B0,TE,theta_B0,R(1));


func=@(b) TotalSignal_Tissue_PerpPlane(sqrt(dchi_a2/b(2)),b(2),B0,TE,theta_B0,R)*b(1)-real(y);
[tmp,t2,t3,flag]=lsqnonlin(func,[rho0,0.4],[0,0.3],[6*rho0,0.6]);

res(1)=tmp(2);
res(2)=sqrt(dchi_a2/tmp(2));

    res(3)=dchi_a2;
    res(4) = tmp(1);
    
%     
%     figure;plot(R,real(y),'ko');
%     hold on;
%     
%     plot(R,func(tmp)+real(y),'r-');
%    disp('');











