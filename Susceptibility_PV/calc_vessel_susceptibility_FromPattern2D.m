function [res,resid,cd_fit,cd,mask]= calc_vessel_susceptibility_FromPattern2D(cd,mask,TE,B0,vox_size,vox_size_interp,lnew,interp,b0,ll,ul)

% lnew: new matrix size[1*2] should be <=size(cd(:,:,1)); for increasing
% % if cd empty: then skip fit, calculate cd_fit with b0.


            FOV = lnew.*vox_size_interp;
            
            nTE=length(TE);
              
            
            cm=roiCOM(mask);
            if ~isempty(cd)
            cd=cd(cm(1)-lnew/2:cm(1)+lnew/2-1,cm(2)-lnew/2:cm(2)+lnew/2-1,:);
            end
            
            mask=mask(cm(1)-lnew/2:cm(1)+lnew/2-1,cm(2)-lnew/2:cm(2)+lnew/2-1);
            
            
                fitfunc2 = @(x) fitfunc(x(1),x(2),x(3),x(4:3+nTE),b0(4+nTE),x(4+nTE:5+nTE),B0,TE,b0(7+nTE),b0(8+nTE),cd,mask,FOV,vox_size,vox_size_interp,interp);
                b0=b0([1:5,7,8]);
                ll=ll([1:5,7,8]);
                ul=ul([1:5,7,8]);

          
          if ~isempty(cd)
             options = optimset('Display','iter','MaxIter',1000);   
             res=b0;
             for i=1:1 
              %res=fminsearch(fitfunc2,res,ll,ul,options);
              tic;
              
              res=fminsearch(fitfunc2,res,options);
            %  res=patternsearch(fitfunc2,res,[],[],[],[],ll,ul,options);
              
            %  res=lsqnonlin(fitfunc2,res,[],[],options);
              fprintf('iter %d: ',i);
              disp(res);
             toc;
             end
          end
             [tmp,cd_fit]=fitfunc2(res);
             if ~isempty(cd)
              resid=sqrt(mean(vec(abs(cd_fit-cd).^2,mask)));  
             else
                 resid=[];
             end
 
 
    
function res=fr0_2D(x,y,rad,dchi,B0,TE,theta,phi_Vessel,Sv,St,bgPhase)
    % theta, phi_Vessel, and bgPhase in degrees;
    
    
    phi=atan2(y,x)*180/pi-phi_Vessel;
    
    r=sqrt(y.^2+x.^2);
   % r=single(r);
    
    phase=dchi2phase(rad,dchi,B0,TE,r,theta,phi)+bgPhase*pi/180;
   
    %phase=s(phase);
    
    S = zeros(size(phase));
    S(r<=rad+0.00001)=Sv;
    S(r>rad+0.00001)=St;
    
    res=S.*exp(1i*phase);
    
function [res,cd_fit]=fitfunc(rad,dchi,Sv,St,bgPhase,center,B0,TE,theta,phi_Vessel,cd,roi,FOV,voxSize,voxSize_interp,interp)
            
res=zeros([sum(roi(:))*2*length(TE),1]);

T2s = 1000/275; %Blockley et al., Magnetic Resonance in Medicine 60:1313–1320 (2008)
for i=1:length(TE)
    
    % when calling fr2; phi is measured relative to the second
    % dimension; now relative to the first dimension.
    
    fr2=@(x,y,z) fr0_2D(x-center(1),y-center(2),rad,dchi,B0,TE(i),theta,phi_Vessel,Sv*exp(-(TE(i)-TE(1))/T2s),St(i),bgPhase);
    
    %cd_fit(:,:,i)=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,true);
    
    cd_fit(:,:,i)=Image_CartesianKSpace_ConvWtFFT(fr2,FOV,voxSize,voxSize_interp,interp);
    
    if ~isempty(cd)
        cd_tmp=cd(:,:,i);
        cd_fit_tmp=cd_fit(:,:,i);
        resid=cd_tmp(roi>0)-cd_fit_tmp(roi>0);       
        res((i-1)*sum(roi(:))*2+1:i*sum(roi(:))*2) = [real(resid(:));imag(resid(:))];
    end
end
           
res=sum(res.^2,1);
           
           
           
           