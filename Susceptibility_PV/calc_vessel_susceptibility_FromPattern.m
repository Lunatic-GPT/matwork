function [res,resid,cd_fit,cd_data,roi]= calc_vessel_susceptibility_FromPattern(cd,roiRad,theta_B0,TE,B0,vox_size,interp_factor,iTE4s,lnew)


% lnew: new matrix size[1*2] should be <=size(cd(:,:,1)).

if length(lnew)==1
    lnew=[lnew,lnew];
end
            roiRad_i=round(mean(roiRad./vox_size.*interp_factor));
            
            
            
            FOV = lnew.*vox_size(1:2)./interp_factor;
            
            nTE=length(iTE4s);
           %fitfunc2 = @(x) fitfunc(x(3),x(4),x(1:2),data(:,:,i),roi3(cropr,cropc),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC);
           %fitfunc(a,dchi,Sv,St,bgPhase,center,B0,TE,theta,cd,roi,FOV,voxSize,voxSize_interp)
           vox_size_new=vox_size./interp_factor;
          for i=1%1:2 % fit a second time with adjusted center position
            if i==1
                center=ceil((size(cd(:,:,1))+1)/2);           
            else
                center=ceil((size(cd(:,:,1))+1)/2)+round(res(end-1:end)./vox_size_new);
            end
            
            cd_data=cd(center(1)-lnew(1)/2:center(1)+lnew(1)/2-1,center(2)-lnew(2)/2:center(2)+lnew(2)/2-1,:);
               
            roi=mask_circle(size(cd(:,:,1)),roiRad_i,center,1);
            roi=roi(center(1)-lnew(1)/2:center(1)+lnew(1)/2-1,center(2)-lnew(2)/2:center(2)+lnew(2)/2-1,:);
            St0=mean_roi(abs(cd_data(:,:,iTE4s)),roi>0)';
          
            if i==1
                b0= [0.1,0.5,St0.*0.5,St0,0,0,0];
                b0=[0.4,0.1,8.9044,1.1604,70,68.2922,0,0,0];
            else
                b0=res;
            end
          
              
          fitfunc2 = @(x) fitfunc(x(1),x(2),x(3:3+nTE-1),x(3+nTE:3+nTE*2-1),x(3+nTE*2),x(4+nTE*2:5+nTE*2),B0,TE(iTE4s),theta_B0,cd_data(:,:,iTE4s),roi,FOV,vox_size,vox_size_new);
            
          % res=lsqnonlin(fitfunc2,[0,0,0.0824,1.3538],[-0.5,-0.5,0,0],[0.5,0.5,1,8]);
      %    b0= [0.1,0.5,St0.*0.5,St0,0,0,0];
          
          ll=[0,0.2,zeros(1,nTE),zeros(1,nTE),-0.1,-0.6,-0.6];
          ul=[1,0.7,St0,St0*2,0.1,0.6,0.6];
          options = optimoptions('lsqnonlin','MaxIteration',1);
          options.Display='iter';
          for iter=1:10
          res=lsqnonlin(fitfunc2,b0,ll,ul,options);
          fprintf('Iter %d: ',iter);
          disp(res);
          b0=res;
          end
          end
                 
          [tmp,cd_fit]=fitfunc2(res);
          resid=sqrt(mean(vec(abs(cd_fit-cd_data(:,:,iTE4s)).^2,roi)));  
          
            
            
function res=fr(r,phi,a,dchi,B0,TE,theta,Sv,St,bgPhase)
    
    phase=dchi2phase(a,dchi,B0,TE,r,theta,phi)+bgPhase;
   
    S = zeros(size(phase));
    S(r<=a)=Sv;
    S(r>a)=St;
    
    res=S.*exp(1i*phase);
    
    
    
    %% 2D 
function [res,cd_fit]=fitfunc(a,dchi,Sv,St,bgPhase,center,B0,TE,theta,cd,roi,FOV,voxSize,voxSize_interp)
            
res=[];

           for i=1:length(TE)
               
               % when calling fr2; phi is measured relative to the second
               % dimension; now relative to the first dimension.
               
            fr2=@(r,phi) fr(r,phi+90,a,dchi,B0,TE(i),theta,Sv(i),St(i),bgPhase);
            
            cd_fit(:,:,i)=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,true);
            
            cd_data=cd(:,:,i);
            cd_fit_tmp=cd_fit(:,:,i);
            resid=cd_data(roi>0)-cd_fit_tmp(roi>0);
            res=cat(1,res,[real(resid(:));imag(resid(:))]);
            
           end
           

           
           
           
           