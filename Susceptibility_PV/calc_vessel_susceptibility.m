function [res,flag,resid_magMoment]= calc_vessel_susceptibility(cd,m_vessel,m_mom,m_background,theta_B0,TE,B0,vox_size,iTE4mom,iTE4s)
% assume B0 projection is along the first dimension
% angles in theta_B0 is in degree
% angle in d is in radian
% TE in ms;
% theta_B0 in degree
% vox_size in mm, one element; assume isotropic;
% res: 1*3; res(1) is dchi (ppm); res(2) is radius in mm; res(3) is dchi*radius^2 (unit: ppm*mm2).

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
mv=ri(m_vessel,'','','d');
m_mom=ri(m_mom,'','','d');
m_background=ri(m_background,'','','d');


m_mom=permute(m_mom,[2,1,3]);
mv=permute(mv,[2,1,3]);
m_background=permute(m_background,[2,1,3]);


m_ind=ind2subb(size(mv),find(mv>0));
center=mean(m_ind,1);

% norm=mean_roi(s4,m2>0&m==0);
% s4=s4/norm(1);
ml=mv;
rho=abs(mean_roi(cd,m_background));
S=mean_roi(cd,ml).*sum(ml(:)>0)*vox_size^2;

R = sqrt(sum(ml(:)>0)/pi)*vox_size;
nTE = length(TE);

if exist('iTE4mom','var') && ~isempty(iTE4mom)
    [dchi_a2,resid_magMoment]=magMoment_outSidePattern(p(:,:,iTE4mom)*180/pi,TE(iTE4mom),B0,theta_B0,m_mom>0,center,vox_size*ones(1,3));
    disp(S);
    if dchi_a2>0
        [res,flag]=vein_QSM_fixMoment(S(iTE4s),rho(iTE4s),TE(iTE4s),R,B0,theta_B0,dchi_a2);
        
        % use the following for debug
        %     [pdata,dchia]=plot_chi_rhoc_vs_signal_fixedm(TE,R,B0,theta_B0,dchi_a2);
        %     %set S = [pick a value on graph]*rho*pi*R*R;
        %     sub=ind2subb(size(pdata),find(abs(real(S/rho/pi/R/R)-real(pdata))<0.0001 & abs(imag(S/rho/pi/R/R)-imag(pdata))<0.0001));
        %      vein_QSM_smallAngle_1TE(S,rho,TE,R,B0,theta_B0,dchi_a2) output
        %      should be the same as dchia(sub(1,1));
        %     hold on;plot(real(S/rho/pi/R/R),imag(S/rho/pi/R/R),'ro');
        res(2) = sqrt(dchi_a2/res);  %in units of mm
        res(3)=dchi_a2;
        p_in = dchi2phase(res(2),res(1),B0,TE,0,theta_B0,0);
    %  res(4) = S/rho;
       res(4:4+nTE-1) = imag(S)./rho/pi/res(2)^2./sin(p_in)';  %normalized rho c.
    
    else
        res=zeros(1,3+nTE);
        flag=-1;
    end
   
   % res(5)=rho(1);
else
    
    [res,flag]=vein_QSM(S(iTE4s),rho(iTE4s),TE(iTE4s),R,B0,theta_B0);
    resid_magMoment=0;
    res(3)=res(1)*res(2)^2;
    p_in = dchi2phase(res(2),res(1),B0,TE,0,theta_B0,0);
    res(4:4+nTE-1) = imag(S)./rho/pi/res(2)^2./sin(p_in)';  %normalized rho c.
  %  res(5)=rho(1);
    
end
fprintf('chi = %3.2f ppm; rad = %3.2f mm; mom = %4.3f ppm.mm2;\n',res(1),res(2),res(3));

fprintf('norm. rhoc = ');
fprintf('%3.2f, ',res(4:end));
fprintf('\n');





