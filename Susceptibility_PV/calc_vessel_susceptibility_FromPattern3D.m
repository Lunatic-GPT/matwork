function fitres = calc_vessel_susceptibility_FromPattern3D(cd0,mask0,TE,B0,vox_size,vox_size_interp,lnew,interp,par)

% lnew: new matrix size[1*2] should be <=size(cd(:,:,1)); for increasing
% if cd empty: then skip fit, calculate cd_fit with b0.
% set ll to [] or not supply if want to use fminsearch, otherwise, use
% patternsearch
% b0 should include [rad,chi,Sv(1),St(1:2),bgphase,theta0,phi0]
global iter;
iter=0;
global fit_vessel_phi;
global fit_method;
global fitphase;
global fit_crossSection;
global skip_fit;
fit_crossSection=true;
fitphase=false;
%fit_method='CG';
%fit_method='fminsearch';
%fit_method = 'twostep';
%fit_method='patternsearch';
fit_method='lsqnonlin';

  FOV = lnew.*vox_size_interp;         
  nTE=length(TE);
              
  vec=thetaPhi2unitVec(par.theta,par.phi);

  if ~fit_crossSection
      nseg=1;
  else
      nseg=floor(par.len/vox_size(1));
  end
      tic;
  for iseg=1:length(par.seg)%1:nseg       
      i=par.seg(iseg);
      fprintf('Seg %d: \n',i);
  
            cm=roiCOM(mask0)+vec*(i-(nseg+1)/2);
            cm=round(cm);
            if ~isempty(cd0)
            cd=cd0(cm(1)-lnew/2:cm(1)+lnew/2-1,cm(2)-lnew/2:cm(2)+lnew/2-1,cm(3)-lnew/2:cm(3)+lnew/2-1,:);
            end
            
            if ~fit_crossSection
              mask=mask0(cm(1)-lnew/2:cm(1)+lnew/2-1,cm(2)-lnew/2:cm(2)+lnew/2-1,cm(3)-lnew/2:cm(3)+lnew/2-1);           
            cm=roiCOM(mask);
            cm=cm-(lnew/2+1)*[1,1,1];
            cm=cm.*vox_size;
            end
 
  if ~fit_crossSection
    b0=[par.rad,par.dchi,cm,par.theta,par.phi];
    ll=[0.2,0.2,-2,-2,-2,par.theta-10,par.phi-10];
    ul=[0.6,0.6,2,2,2,par.theta+10,par.phi+10];
    fitfunc2=@(x) fitfunc(x,par,B0,TE,cd,mask,FOV,vox_size,vox_size_interp,interp);
       if i>1
           break;
       end
  else
    b0=[par.rad,par.dchi,par.cm];
    ll=[0.1,0.01,-par.roi_rad,-par.roi_rad];
    ul=[1,0.8,par.roi_rad,par.roi_rad];
    fitfunc2=@(x) fitfunc_rot(x,par,B0,TE,cd,FOV,vox_size,vox_size_interp,interp,par.roi_rad);   
    if fit_vessel_phi
       b0(end+1)=par.phi;
       ll(end+1)=par.phi-20;
       ul(end+1)=par.phi+20;    
    end
    
    if par.fit_S
       b0(end+1:end+3)=[par.Sv,par.St];
       ll(end+1:end+3)=[par.Sv,par.St]*0.5;
       ul(end+1:end+3)=[par.Sv,par.St]*1.5;    
    end
    
  end
  
  
          if isempty(cd)
              
              [tmp,cd_fit,cd_rot,cd_fit_rot,roi_circ]=fitfunc2(b0);
                  fitres.cd_fit_rot=cd_fit_rot;
                  fitres.cd_rot=cd_rot;
                  fitres.roi_circ=roi_circ;
                  fitres.cd_fit=cd_fit;
                  
          else
              %res=fminsearch(fitfunc2,res,ll,ul,options);
            
              if strcmp(fit_method,'lsqnonlin')
                  rls=version('-release');
                  rls=str2num(rls(1:4));
                 
                      if rls<=2013
                          options = optimoptions('lsqnonlin','Display','iter','MaxIter',Inf,'TolFun',0,'TolPCG',0,'MaxFunEvals',Inf,'FinDiffRelStep',0.02);
                      else
                          options = optimoptions('lsqnonlin','Display','iter','MaxIter',Inf,...
                              'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,'OptimalityTolerance',0,'StepTolerance',1e-6,...
                              'FiniteDifferenceStepSize',2e-2);
                      end
                      
                      if ~skip_fit
                          [res,RESNORM,resid,exitflag,output]=lsqnonlin(fitfunc2,b0,ll,ul,options);
                          
                 
                      else
                          res=b0;
                          resid=0;
                          exitflag=0;
                          output=0;
                         % fitfunc2(b0);
                      end
                  
              elseif strcmp(fit_method,'CG')
                  options = optimoptions('lsqnonlin','Display','iter','MaxIter',Inf,'StepTolerance',1e-6,'FiniteDifferenceStepSize',1e-2);
                  res= CG_backTrack(fitfunc2,b0,options);
%               elseif strcmp(fit_method,
%                   
              elseif strcmp(fit_method,'fminsearch')
                  
               options = optimset('Display','iter','MaxIter',Inf);
               [res,resid,exitflag,output]=fminsearch(fitfunc2,b0,options);
             
              else
                options = optimset('Display','iter','MaxIter',1000);
                [res,resid,exitflag,output]=patternsearch(fitfunc2,b0,[],[],[],[],ll,ul,options);
              end
              if  ~fit_crossSection
                  [tmp,cd_fit]=fitfunc2(res);
                  
              else
                  [tmp,cd_fit,cd_rot,cd_fit_rot,roi_circ]=fitfunc2(res);
                  fitres(iseg).cd_fit_rot=cd_fit_rot;
                  fitres(iseg).cd_rot=cd_rot;
                  fitres(iseg).roi_circ=roi_circ;
              end
              
              fitres(iseg).res=res;
              fitres(iseg).resid=resid;
              fitres(iseg).cd_fit=cd_fit;
              fitres(iseg).cd=cd;            
              fitres(iseg).exitflag=exitflag;   
              fitres(iseg).output=output;
                    disp(res);
              %% for debug
              debug=0;
              if debug
                  if ~fit_crossSection
                     
                      if ~fitphase
                          showfit(cd,cd_fit,mask);
                      else
                          showfit(angle(cd),angle(cd_fit),mask);
                      end
                      cd_fit_rot=[];
                      cd_rot=[];
                  else
                    
                      if fitphase
                          showfit(angle(cd_rot(:,:,floor(end/2)+1,:)),angle(cd_fit_rot(:,:,floor(end/2)+1,:)),roi_circ);
                          
                      else
                          showfit(cd_rot(:,:,floor(end/2)+1,:),cd_fit_rot(:,:,floor(end/2)+1,:),roi_circ);
                          
                      end
                  end
                  
                
              end
          end
        fprintf('Elapsed/remaining times: %5.1f/%5.1f\n',toc/60,toc/60*(nseg-i)/i);  
  end  
    function showfit(cd,cd_fit,roi)
        
        
        global fitphase;  
        ind=roi_ind2sub(roi);
        if size(roi,3)>1
            sl=unique(ind(:,3));
        else
            sl=1;
        end
        cd=setv_roi(cd,roi==0,0);
        cd_fit=setv_roi(cd_fit,roi==0,0);
        
        if fitphase
            figure(100);imshow4(cat(4,cd(:,:,sl,:),cd_fit(:,:,sl,:)),[-pi/3,pi/3],[2*size(cd,4),length(sl)]);
            
            
            figure(102);plot(vec(cd_fit,roi),vec(cd,roi),'o');
            hold on;plot(min_max(vec(cd_fit,roi)),min_max(vec(cd_fit,roi)),'r-');
        else
            figure(100);imshow4(cat(3,real(cd(:,:,sl)),real(cd_fit(:,:,sl))),[],[2,length(sl)]);
            
            figure(101);imshow4(cat(3,imag(cd(:,:,sl)),imag(cd_fit(:,:,sl))),[],[2,length(sl)]);
            
            figure(102);plotc(vec(cd,roi),vec(cd_fit,roi),'o');
            
            xx=[vec(real(cd),roi);vec(imag(cd),roi);vec(real(cd_fit),roi);vec(imag(cd_fit),roi)];
           % xx=vec([real(cd_fit(:)),imag(cd_fit(:))],roi);
            hold on;plot(min_max(xx),min_max(xx),'r-');

        end
        
function [res,cd_fit]=fitfunc(x,par,B0,TE,cd,roi,FOV,voxSize,voxSize_interp,interp)


    nTE = length(TE);
    Sv=par.Sv;
    St=par.St;
    dchi=x(2);
    rad=x(1);
    
    
    bgPhase=par.bgPhase;
    center=x(3:5);
    theta=x(6);
    phi_Vessel=x(7);
    
    T2s_v=par.T2s_v;
 

global fitphase;

if fitphase
 res=zeros([sum(roi(:))*nTE,1]);
else
 res=zeros([sum(roi(:))*2*nTE,1]);
end

for i=1:nTE
    
    % when calling fr2; phi is measured relative to the second
    % dimension; now relative to the first dimension.
    
    fr2=@(x,y,z) fr(x-center(1),y-center(2),z-center(3),rad,dchi,B0,TE(i),theta,phi_Vessel,Sv*exp(-(TE(i)-TE(1))/T2s_v),St(i),bgPhase);
    
    %cd_fit(:,:,i)=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,true);
   
    cd_fit(:,:,:,i)=Image_CartesianKSpace_ConvWtFFT(fr2,FOV,voxSize,voxSize_interp,interp);
   
    if ~isempty(cd)
        cd_tmp=cd(:,:,:,i);
        cd_fit_tmp=cd_fit(:,:,:,i);
        if fitphase       
           ph_resid=conj(cd_tmp).*cd_fit_tmp; 
           res((i-1)*sum(roi(:))+1:i*sum(roi(:))) = angle(ph_resid(roi>0));
        else
           resid=cd_tmp(roi>0)-cd_fit_tmp(roi>0);
           res((i-1)*sum(roi(:))*2+1:i*sum(roi(:))*2) = [real(resid(:));imag(resid(:))];
        end
    else
        res=0;
    end
    
end

global fit_method;

if ~strcmp(fit_method,'lsqnonlin')

res=sum(res.^2,1);
end
           
function [res,cd_fit,cd_rot,cd_fit_rot,roi_circ]=fitfunc_rot(x,par,B0,TE,cd,FOV,voxSize,voxSize_interp,interp,roi_rad)

tic;
    nTE = length(TE);
   

    if isfield(par,'fixchi') && par.fixchi
       dchi=par.dchi;
       rad=sqrt(x(1)/100/dchi); 
       np=1;
    elseif isfield(par,'fixmom') && par.fixmom
       dchi=x(1);
       np=1;
       rad=sqrt(par.dchi*par.rad^2/dchi);
    else    
     dchi=x(2);
     rad=x(1);
     np=2;
    end
    
    bgPhase=par.bgPhase;
    
global fit_vessel_phi;  
    theta=par.theta;
     
     if fit_vessel_phi
        phi_Vessel = x(np+1);
        np=np+1;
     else
        phi_Vessel=par.phi;
     end
     
    mat=transform_matrix_rotation(theta,phi_Vessel);
    global fit_crossSection;
    if ~fit_crossSection
        center=x(np+1:np+3);
        np=np+3;
    else
        center=mat*[x(np+1);x(np+2);0];%x(3) moving up ; x(4) moving left.
        np=np+2;
    end
    
    %T2s_v=par.T2s_v;
     T2s_v=dchi2T2(dchi*1e-6)*1000;
     if ~par.fit_S
         
       
       St=par.St;
    %Sv=par.Sv;
       Sv= par.Swm_Ref(1).*exp(-TE(1)/T2s_v).*exp(TE(1)/par.T2s_t)*1.05;
       
     else
        Sv=x(np+1);
        St=x(np+2:np+3);
        np=np+3;
     end
    
     
global fitphase;

res=[];

cm=floor(size(cd(:,:,:,1))/2)+1;

vec=thetaPhi2unitVec(theta,par.phi);
p1=cm+vec;
p2=cm-vec;

%mask_circle(dim,rad,center,include_equal)
for i=1:nTE
    
    % when calling fr2; phi is measured relative to the second
    % dimension; now relative to the first dimension.
    
    fr2=@(x,y,z) fr(x-center(1),y-center(2),z-center(3),rad,dchi,B0,TE(i),theta,phi_Vessel,Sv*exp(-(TE(i)-TE(1))/T2s_v),St(i),bgPhase);
      
    %cd_fit(:,:,i)=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,true);
    cd_fit(:,:,:,i)=Image_CartesianKSpace_ConvWtFFT(fr2,FOV,voxSize,voxSize_interp,interp);  
    B0_dir=[0,0,1];
    crop_size = size(cd_fit(:,:,:,i));
    
      [mg,p,deg,tmp,fov_center]=rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,voxSize_interp,cd_fit(:,:,:,i),crop_size,1,voxSize_interp);
    %  roi_center=ceil((size(tmp)+1)/2);
      cd_fit_rot(:,:,i)=tmp(:,:,fov_center(3));
    
      
      
    if ~isempty(cd)
        
      
      [mg,p,deg,tmp]=rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,voxSize_interp,cd(:,:,:,i),crop_size,1,voxSize_interp);

      
      cd_rot(:,:,i)=tmp(:,:,fov_center(3));
      
      
      roi_rad_i=round(roi_rad/voxSize_interp(1));
        %cm_new=centerFinder(cd_rot,1,roi_rad_i,1);
      cm_new=fov_center;
      roi_circ=mask_circle(size(cd_rot(:,:,1,1)),roi_rad_i,cm_new(1:2),1);

        cd_tmp=cd_rot(:,:,i);
        cd_fit_tmp=cd_fit_rot(:,:,i);
        if fitphase       
           ph_resid=conj(cd_tmp(roi_circ>0)).*cd_fit_tmp(roi_circ>0); 
           res=cat(1,res,angle(ph_resid(:)));
        else
          resid=cd_tmp(roi_circ>0)-cd_fit_tmp(roi_circ>0);
       %   res((i-1)*sum(roi(:))*2+1:i*sum(roi(:))*2) = [real(resid(:));imag(resid(:))];         
          res=cat(1,res,real(resid(:)),imag(resid(:)));
         end
    else
        res=0;
    end
    
end
res2=sum(res.^2,1);

%res=double(res);
global fit_method;

% if strcmp(fit_method,'lsqnonlin')
%     global iter;
%     
%     if mod(iter,5)==0
%         save tmp x res2
%     else
%         i=mod(iter,5);
%         tmp=load('tmp.mat');
%         grad=(tmp.res2-res2)/(tmp.x(i)-x(i));
%         fprintf('%8.6f   ',[i,grad]);
%         
%     end
%     iter=iter+1;
% else
   fprintf('%8.6f   ',[toc,x,res2,rad,dchi]);
% end
fprintf('\n');

if ~strcmp(fit_method,'lsqnonlin')

res=sum(res.^2,1);
end             

function res=fr(x,y,z,rad,dchi,B0,TE,theta,phi_Vessel,Sv,St,bgPhase)
    % theta, phi_Vessel, and bgPhase in degrees;
    
   % phi = acos(z./sqrt(x.^2+y.^2+z.^2));
    
    
    unitVec=thetaPhi2unitVec(theta,phi_Vessel);
    
    lin=cat(1,-0.5*unitVec,0.5*unitVec);
    
    xyz=cat(4,x,y,z);
    xyzp=Projection2Plane(xyz,unitVec);
    r = distance2Line(lin,xyz);

    clear xyz;
    
    
    zp=Projection2Plane([0,0,1],unitVec);
    phi = angle_bw_2vec(xyzp,zp);
    clear xyzp;
    
    phi(isnan(phi))=0;
    
    %r=single(r);
    
    phase=dchi2phase(rad,dchi,B0,TE,r,theta,phi)+bgPhase*pi/180;
    clear phi;
  %  phase=single(phase);
  
    S = zeros(size(phase));
    S(r<=rad+0.00001)=Sv;
    S(r>rad+0.00001)=St;
      
%     S = zeros(size(phase),'single');
%     S(r<=rad)=Sv;
%     S(r>rad)=St;
%     
    res=S.*exp(1i*phase);
 
    % fitfunc(x,par,B0,TE,cd,mask,FOV,vox_size,vox_size_interp,interp);