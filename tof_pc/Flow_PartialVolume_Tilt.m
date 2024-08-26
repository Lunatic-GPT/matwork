function Flow_PartialVolume_Tilt(params,save_name)
% params should contain the following fields:
% roi: vessel mask
% m_wm: white matter mask
% voxSize: acquired voxSize; (cm)
% voxSize_interp: voxSize of the images; (cm)
% thk: slice thickness: (cm)
% PhaseDifference: whether the phase image is for phase difference between
%             on and off, or for the raw phases.  In the latter case, there
%             will be nroi*2 phase images since the phase image will be 
%             calculated separately for each roi.
% mag: magnitude images;
% ph: phase images; (value should be in degrees).
% FA: flip angle;
% roiRad: radius for the circular ROI encompassing all pixels to be included in model image fitting (cm).
% bgRad: background ROI out-radius for calculating tissue signal (cm)
% interp: interpolation factor for PSF convolution.
% VENC: (cm/s)
% T1: T1 for blood and WM tissue, respectively. (s)
% T2: T2 for blood and WM tissue, respectively. (s)
% TR: (s)
% TE: (s)
% negPhase: whether the flow produce a negative phase
% phi_tilt and theta_tilt: (degree) polar angles for vessel orientation.
% %                         Dimension: 1*nroi. 
%                           x,y,z are the 1st,2nd and 3rd dimensions,
%                           respectively. 
%                           Axes point to the directions of increasing
%                           index value for x and y, and along the flow for
%                           z
%                            
%                           theta_tilt: angle between vessel (pointing toward flow) and +z;
%                           [0,90].
%                           phi_tilt: angle between vessel and +x in the
%                           x-y plane.
% plot_fit: plot fit results
% skip_fit: skip fitting
% v0, rad0: radius and flow velocity in mm and cm/s; only used when
% skip_fit is true;

if ~isfield(params,'plot_fit')
    plot_fit=false;
else
    plot_fit=params.plot_fit;
end

if ~isfield(params,'skip_fit')
    skip_fit=false;
else
    skip_fit=params.skip_fit;
end

roi=params.roi;

mag=params.mag;
ph=params.ph;
% FOV=par.FOV;
% voxSize=par.voxSize;
% voxSize_interp=par.voxSize_interp;
% interp=par.interp;
% VENC=par.VENC;

% bgRad
% num2deg
%params.do_detrend
%params.neg_phase
% params.PhaseDifference
% if true, the phase image should only have one volume,
% if false: the phase image should have nroi*2 volumes;  two for each roi,
% where nroi=max(roi(:));

if max(roi(:))==1  %not yet clusterized
    roi=clusterize2(roi);
end

nroi=max(roi(:));
pv=zeros(1,nroi);
v=zeros(1,nroi);
rad = zeros(1,nroi);
nv=zeros(1,nroi);
ph_total=zeros(1,nroi);
ph_max=zeros(1,nroi);

center=zeros(2,nroi);

pix_area=prod(params.voxSize_interp);



if params.negPhase
    ph=-ph;
end


for i=1:max(roi(:))
    %%
    
    
    tic;
    
    if isfield(params,'vessels') && ~any(i==params.vessels)
        continue;
    end
    
    nv(i)=sum(roi(:)==i);
    if nv(i)==0 || i> length(params.phi_tilt) || isnan(params.phi_tilt(i)) % skip if orientation not determined or no vessel found
        continue;
    end
 %data = preprocess_data(mag,ph,roi,params)
    if ~params.PhaseDifference
        cd_data(:,:,:,i)=preprocess_data(mag(:,:,:,2*i-1:2*i),ph(:,:,:,2*i-1:2*i),roi==i,params);
    else
        cd_data(:,:,:,i)=preprocess_data(mag,ph,roi==i,params);
    end
    
    roiRad_i = (size(cd_data,1)-1)/2;
    roi2=mask_circle([2*roiRad_i+1,2*roiRad_i+1],roiRad_i,[roiRad_i+1,roiRad_i+1],1);
    
  [za,va,sa,st] = prep_interp(params,i);
  
    par_fit=params;
    par_fit.za=za;
    par_fit.sa=sa;
    par_fit.va=va;
    par_fit.st=st;        
    par_fit.roi=roi2;
    par_fit.phi_tilt=params.phi_tilt(i);
    par_fit.theta_tilt=params.theta_tilt(i);
    
    
    fitfunc2 = @(x) fitfunc(x,par_fit,cd_data(:,:,:,i)); 
    %       options = optimoptions('patternsearch','Display','iter','MaxIter',Inf,'MaxFunctionEvaluations',Inf);
    %      [res,tmp,exitflag(i),output(i)]=patternsearch(fitfunc2,[0,0,0.08,0.6],[],[],[],[],[-0.5,-0.5,0.01,0],[0.5,0.5,0.25,5],[],options);
    if datenum(version('-date'))>=736580
    options = optimoptions('lsqnonlin','MaxIter',Inf,'StepTolerance',1e-4,...
        'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,'OptimalityTolerance',0,'Display','iter','FiniteDifferenceStepSize',0.01);
    else
        options = optimoptions('lsqnonlin','MaxIter',Inf,'TolX',1e-4,...
        'MaxFunEvals',Inf,'TolFun',0,'TolPCG',0,'Display','iter','FinDiffRelStep',0.01);
    end

% center in cm
% and radius in mm 
% velocity in 10 cm/s
% phase in 100 deg

    ll=[-0.04,-0.04,0.01,0,(par_fit.phi_tilt-90)/100];
    ul=[0.04,0.04,0.25,0.399,(par_fit.phi_tilt+90)/100];
    
    if skip_fit
        
        res=[0,0,params.rad0,params.v0/10,params.phi_tilt/100];
    else
    [res,tmp,tmp2,exitflag(i),output(i)]=lsqnonlin(fitfunc2,[0,0,0.1,0.1,par_fit.phi_tilt/100],ll,ul,options);
    end
    
    pv(i) = pi*res(3)^2/nv(i)/pix_area;
    v(i) = res(4);  % in 10 cm/s
    rad(i)=res(3); %% in mm
    center(:,i)=res(1:2);
    
    ph_total(i)=mean(ph(roi==i));
    ph_max(i)=max(ph(roi==i));
    
    [resid,cd_fit(:,:,:,i)]=fitfunc2(res);
    
  
    if plot_fit  %check fit
        %disp(res);
        tmp1=cd_data(:,:,i);
        tmp1(~roi2)=0;
        tmp2=cd_fit(:,:,i);
        tmp2(~roi2)=0;
        
        figure; compare_cd(tmp1,tmp2);
        
        fprintf('radius = %3.2f mm; v = %3.2f cm/s\n',res(3),res(4)*10);
    end
    
    fprintf('Time remaining %4.3f s\n',toc*(sum(~isnan(params.theta_tilt))-i));
end
% 
% flow=pv.*v*pix_area.*nv*1000;  %in mm3/s
% flow_uc=ph_total*VENC/180*pix_area.*nv*1000;
% flow_max=ph_max*VENC/180*pix_area*1000;


params=rmfield(params,{'mag','ph','roi','m_wm'});
save(save_name,'v','rad','cd_fit','cd_data','center','params','exitflag','output');
  

function [za,va,sa,st]=prep_interp(params,i)

% i is the vessel roi value;

    theta_tilt=params.theta_tilt(i);

    f_prep=param2name('FPV',theta_tilt,params.TR,params.thk,params.FA,params.T1,params.TE,params.T2);
    
    d_prep=fileparts(mfilename('fullpath'));
    d_prep=fullfile(d_prep,f_prep);
    
    if ~exist(d_prep,'file')
        
        za = linspace(-0.3,0.3,400); %x
        va=linspace(0,8,400); % y
        sa=zeros(length(va),length(za));
        st=zeros(1,length(za));
        for j=1:length(va)
            
            for k=1:length(za)              
                sa(j,k)=signal_flow(za(k),theta_tilt,va(j),params.TR,params.thk,params.FA,params.T1(1))*exp(-params.TE/params.T2(1));
            end
        end
        
        for k=1:length(za)
           st(k) = signal_flow(za(k),0,0,params.TR,params.thk,params.FA,params.T1(2))*exp(-params.TE/params.T2(2));%same as the steady state signal
        end
        save(d_prep,'za','va','sa','st');
        
    else
        load(d_prep);
    end
  
function data = preprocess_data(mag,ph,roi,params)

mag=double(mag);
voxSize_interp=params.voxSize_interp;
bgRad=params.bgRad;
roiRad=params.roiRad;
m_wm=params.m_wm;
PhaseDifference=params.PhaseDifference;

    bgRad_i(1)=round(mean(bgRad(1)./voxSize_interp));
    bgRad_i(2)=round(mean(bgRad(2)./voxSize_interp));
    
    roiRad_i=round(mean(roiRad./voxSize_interp));  % 1 and 2 should be equal
    
pos=roiCOM(roi);

bg_out=mask_circle(size(roi),bgRad_i(1),pos,1);
bg_inn=mask_circle(size(roi),bgRad_i(2),pos,1);

bgroi=bg_out&bg_inn==0&m_wm>0;

mbg=mean_roi(mag,bgroi);

cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;
cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;


ph=tof_bg_rm(ph,bgroi,[],[],false);
mag=tof_bg_rm(mag,bgroi);
        
  

mag1=mag(:,:,:,1);
mag2=mag(:,:,:,2);

if PhaseDifference
    cd=mag2.*exp(1i*(ph)/180*pi) - mag1;% negative phase
else
    cd=mag2.*exp(-1i*(ph(:,:,:,2))/180*pi) - mag1.*exp(-1i*(ph(:,:,:,1))/180*pi);% negative phase
end

cd=cd/mbg(1);
data=cd(cropr,cropc);




function [res,cd_fit,im0,im1]=fitfunc(x,par,cd_data)
% step size 0.001
% center in cm
% and radius in mm 
% velocity in 10 cm/s
% phase in 100 deg

t_fit=tic;
voxSize=par.voxSize;
voxSize_interp=par.voxSize_interp;

FOV=size(cd_data).*voxSize_interp;
interp=par.interp;  % interpolation factor for convolution.
% VENC=par.VENC;
% theta_tilt=par.theta_tilt;
%phi_tilt=par.phi_tilt;
roi=par.roi;
za=par.za;
va=par.va;
sa=par.sa;
st=par.st;

center=x(1:2);
rad=x(3)*0.1;
vmean=x(4)*10;
phi_tilt=x(5)*100;


disp([center,rad,vmean,phi_tilt]);

% 
% % tilt in degree
% if tilt>0

st_int=sum(st)*(za(2)-za(1));
ftissue=@(x,y) st_int*ones(size(x));
imt=Image_CartesianKSpace_ConvWtFFT(ftissue,FOV,voxSize,voxSize_interp,interp); % all pixel value same as ftissue;

par2=par;
par2.VENC=Inf;



if ~par.consider_shift_due_to_PE
    ftrue1=@(x,y) integrate_signal_flow2(x-center(1),y-center(2),vmean,za,va,sa,rad,par);
    ftrue0=@(x,y) integrate_signal_flow2(x-center(1),y-center(2),vmean,za,va,sa,rad,par2);
     
        
im1=Image_CartesianKSpace_ConvWtFFT(ftrue1,FOV,voxSize,voxSize_interp,interp);
im0=Image_CartesianKSpace_ConvWtFFT(ftrue0,FOV,voxSize,voxSize_interp,interp);
else
    
    zi=linspace(-0.3,0.3,100);
     
    for i=1:length(zi)
 
    ftrue1=@(x,y) signal_flow2(x-center(1),y-center(2),zi(i),vmean,za,va,sa,rad,par);  % do it with interp
    ftrue0=@(x,y) signal_flow2(x-center(1),y-center(2),zi(i),vmean,za,va,sa,rad,par2);  % do it with interp
    fshift=@(x,y) shift_flow(x-center(1),y-center(2),zi(i),vmean,rad,par);

    im1(:,:,i)=Image_CartesianKSpace_ConvWtFFT_FlowShiftPE(ftrue1,fshift,FOV,voxSize,voxSize_interp,interp);
    im0(:,:,i)=Image_CartesianKSpace_ConvWtFFT_FlowShiftPE(ftrue0,fshift,FOV,voxSize,voxSize_interp,interp);
    
    
    end
    
    im1=sum(im1,3)*(zi(2)-zi(1));
    im0=sum(im0,3)*(zi(2)-zi(1));
    
    
    
end
cd_fit=(im1-im0)*par.lambda/mean_roi(imt,roi>0);%

if exist('cd_data','var') && ~isempty(cd_data)
    resid=cd_data(roi>0)-cd_fit(roi>0);
    res=[real(resid(:));imag(resid(:))];
else
    res=0;
end


toc(t_fit);
% global     use_fminsearch;
% if use_fminsearch
%     res=sum(res.^2,1);
% end

function res=shift_flow(x,y,z,vmean,rad,params)  % do it with interp

if numel(z)==1
    z=repmat(z,size(x));
end

theta_tilt=params.theta_tilt;
phi_tilt=params.phi_tilt;


r2=xyz2r(x,y,z,theta_tilt,phi_tilt);

sel=r2<=rad*rad;
vi=vmean*2*(1-r2(sel)/rad/rad);  %

res=r2*0;
res(sel)=vi*sin(theta_tilt*pi/180)*sin(phi_tilt*pi/180)*params.Time2TE;






function res=integrate_signal_flow2(x,y,vmean,za,va,sa,rad,params)

res=0*x;
phi_tilt=params.phi_tilt;
cphi=cos(phi_tilt*pi/180);
sphi=sin(phi_tilt*pi/180);
%stheta=tan(theta_tilt*pi/180);
%ctheta=cos(theta_tilt*pi/180);
%debug=0*x;
%r2=repmat(0*x,[1,1,50]);

%ztmp=linspace(-0.3,0.3,50);
for i=1:size(x,1)

   
    for j=1:size(x,2)
        
         v=[x(i,j),y(i,j)];
    %     vtilt=[cphi,sphi];
         vperp=[sphi,-cphi];
         
         if abs(sum(v.*vperp,2))>rad
             continue;
         end
         
     %   debug(i,j)=1;
       % tic;
        f=@(z) signal_flow2(x(i,j),y(i,j),z,vmean,za,va,sa,rad,params);
      % disp([i,j]);
      
     % tic;
     if 0 % alternative method
        res(i,j)=integral(f,-0.3,0.3,'RelTol',1e-3,'AbsTol',1e-6);
      %  toc;
     else
       % tic;
        zi=linspace(-0.3,0.3,100);
        tmp=f(zi);
        res(i,j)=mean(tmp)*0.6;
     end
        %toc;
       % zi=linspace(-0.3,0.3,1000);
       % tmp=f(
        %res(i,j)=
       %toc;
   %   [tmp,r2(i,j,:)]=signal_flow2(x(i,j),y(i,j),ztmp,vmean,za,va,sa,rad,params);
      
    end
   
end

disp('');

    %     f=@(z) signal_flow2(x(:),y(:),z,vmean,za,va,sa,rad,VENC,theta_tilt,phi_tilt);
    %     res=integral(f,-0.3,0.3,'ArrayValued',true);

function [res,r2]=signal_flow2(x,y,z,vmean,za,va,sa,rad,params)  % do it with interp

if numel(z)==1
    z=repmat(z,size(x));
end

VENC = params.VENC;
theta_tilt=params.theta_tilt;
phi_tilt=params.phi_tilt;


r2=xyz2r(x,y,z,theta_tilt,phi_tilt);

sel=r2<=rad*rad;
vi=vmean*2*(1-r2(sel)/rad/rad);  %


%  tint=tic;
z2=z(sel)-vi*params.TE*cos(theta_tilt*pi/180);
z2(z2<za(1))=za(1);
res_tmp= interp2(za,va,sa,z2,vi);
res_tmp=res_tmp.*exp(1i*vi/VENC*pi*cos(theta_tilt*pi/180));

res=0*r2;
res(sel)=res_tmp;

function d2=xyz2r(x,y,z,theta_tilt,phi_tilt)
% tilt in degree;
theta_tilt=theta_tilt*pi/180;
phi_tilt=phi_tilt*pi/180;


lin=[0,0,0;sin(theta_tilt)*cos(phi_tilt),sin(theta_tilt)*sin(phi_tilt),cos(theta_tilt)];

%z=z*ones(size(x));

d2=distance_sq_2Line(lin,x,y,z);


function res=distance_sq_2Line(lin,x,y,z)



c1=x-lin(1,1);
c2=y-lin(1,2);
c3=z-lin(1,3);

b=lin(2,:)-lin(1,:);
%ac=sos(c,2);
ac=sqrt(c1.^2+c2.^2+c3.^2);
ab=sos(b,2);



res=(1-((c1*b(1)+c2*b(2)+c3*b(3))./ac/ab).^2).*ac.*ac;


%%

res(isnan(res))=0;



function sart=signal_flow(z,tilt,v,TR,thk,fa,T1)
% pos and fa gives the RF slice selection profile;
% for an ideal boxcar profile use pos=[-thk/2 thk/2], fa=[fa,fa];
% assuming the vessel center pass through [0,0,0];
% tilt angle in degree from norm toward x axis in the x-z plane;
% right-hand coordinate system.
% tilt in degree
% v in cm/s
% thk: in cm
% Time in s.
if numel(v)==1
    v=v*ones(size(z));
end

tilt=tilt/180*pi;

sz_init=1;
pos= linspace(-6,6,400);
faNorm=[3.0e-05 3.0e-05 2.9e-05 2.5e-05 2.0e-05 1.4e-05 5.9e-06 2.3e-06 1.1e-05 1.9e-05 2.6e-05 3.2e-05 3.6e-05 3.8e-05 3.8e-05 3.5e-05 3.0e-05 2.3e-05 1.5e-05 4.8e-06 5.7e-06 1.6e-05 2.6e-05 3.5e-05 4.2e-05 4.7e-05 4.9e-05 4.8e-05 4.4e-05 3.7e-05 2.7e-05 1.6e-05 2.6e-06 1.1e-05 2.5e-05 3.7e-05 4.8e-05 5.7e-05 6.2e-05 6.3e-05 6.1e-05 5.5e-05 4.5e-05 3.2e-05 1.6e-05 1.7e-06 2.0e-05 3.8e-05 5.4e-05 6.8e-05 7.8e-05 8.4e-05 8.5e-05 8.0e-05 7.1e-05 5.6e-05 3.7e-05 1.5e-05 9.8e-06 3.5e-05 5.9e-05 8.1e-05 9.9e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.3e-05 7.1e-05 4.3e-05 1.0e-05 2.5e-05 6.1e-05 9.5e-05 1.2e-04 1.5e-04 1.6e-04 1.7e-04 1.7e-04 1.5e-04 1.3e-04 9.0e-05 4.6e-05 4.3e-06 5.8e-05 1.1e-04 1.6e-04 2.1e-04 2.4e-04 2.6e-04 2.7e-04 2.6e-04 2.3e-04 1.9e-04 1.3e-04 5.5e-05 2.6e-05 1.1e-04 2.0e-04 2.8e-04 3.4e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.3e-04 2.4e-04 1.2e-04 2.3e-05 1.8e-04 3.5e-04 5.1e-04 6.6e-04 7.9e-04 8.8e-04 9.3e-04 9.3e-04 8.7e-04 7.5e-04 5.7e-04 3.3e-04 5.0e-05 2.7e-04 6.0e-04 9.2e-04 1.2e-03 1.4e-03 1.6e-03 1.6e-03 1.5e-03 1.3e-03 8.7e-04 2.7e-04 4.9e-04 1.4e-03 2.4e-03 3.6e-03 4.7e-03 5.7e-03 6.6e-03 7.2e-03 7.5e-03 7.1e-03 6.1e-03 4.1e-03 1.2e-03 3.0e-03 8.6e-03 1.6e-02 2.5e-02 3.5e-02 4.8e-02 6.3e-02 8.0e-02 1.0e-01 1.2e-01 1.5e-01 1.7e-01 2.0e-01 2.3e-01 2.6e-01 3.0e-01 3.3e-01 3.7e-01 4.1e-01 4.4e-01 4.8e-01 5.2e-01 5.6e-01 6.0e-01 6.4e-01 6.7e-01 7.1e-01 7.4e-01 7.7e-01 8.0e-01 8.3e-01 8.6e-01 8.8e-01 9.0e-01 9.2e-01 9.4e-01 9.5e-01 9.7e-01 9.8e-01 9.8e-01 9.9e-01 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 9.9e-01 9.8e-01 9.8e-01 9.7e-01 9.5e-01 9.4e-01 9.2e-01 9.0e-01 8.8e-01 8.6e-01 8.3e-01 8.0e-01 7.7e-01 7.4e-01 7.1e-01 6.7e-01 6.4e-01 6.0e-01 5.6e-01 5.2e-01 4.8e-01 4.4e-01 4.1e-01 3.7e-01 3.3e-01 3.0e-01 2.6e-01 2.3e-01 2.0e-01 1.7e-01 1.5e-01 1.2e-01 1.0e-01 8.0e-02 6.3e-02 4.8e-02 3.5e-02 2.5e-02 1.6e-02 8.6e-03 3.0e-03 1.2e-03 4.1e-03 6.1e-03 7.1e-03 7.5e-03 7.2e-03 6.6e-03 5.7e-03 4.7e-03 3.6e-03 2.4e-03 1.4e-03 4.9e-04 2.7e-04 8.7e-04 1.3e-03 1.5e-03 1.6e-03 1.6e-03 1.4e-03 1.2e-03 9.2e-04 6.0e-04 2.7e-04 5.0e-05 3.3e-04 5.7e-04 7.5e-04 8.7e-04 9.3e-04 9.3e-04 8.8e-04 7.9e-04 6.6e-04 5.1e-04 3.5e-04 1.8e-04 2.3e-05 1.2e-04 2.4e-04 3.3e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.4e-04 2.8e-04 2.0e-04 1.1e-04 2.6e-05 5.5e-05 1.3e-04 1.9e-04 2.3e-04 2.6e-04 2.7e-04 2.6e-04 2.4e-04 2.1e-04 1.6e-04 1.1e-04 5.8e-05 4.3e-06 4.6e-05 9.0e-05 1.3e-04 1.5e-04 1.7e-04 1.7e-04 1.6e-04 1.5e-04 1.2e-04 9.5e-05 6.1e-05 2.5e-05 1.0e-05 4.3e-05 7.1e-05 9.3e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.9e-05 8.1e-05 5.9e-05 3.5e-05 9.8e-06 1.5e-05 3.7e-05 5.6e-05 7.1e-05 8.0e-05 8.5e-05 8.4e-05 7.8e-05 6.8e-05 5.4e-05 3.8e-05 2.0e-05 1.7e-06 1.6e-05 3.2e-05 4.5e-05 5.5e-05 6.1e-05 6.3e-05 6.2e-05 5.7e-05 4.8e-05 3.7e-05 2.5e-05 1.1e-05 2.6e-06 1.6e-05 2.7e-05 3.7e-05 4.4e-05 4.8e-05 4.9e-05 4.7e-05 4.2e-05 3.5e-05 2.6e-05 1.6e-05 5.7e-06 4.8e-06 1.5e-05 2.3e-05 3.0e-05 3.5e-05 3.8e-05 3.8e-05 3.6e-05 3.2e-05 2.6e-05 1.9e-05 1.1e-05 2.3e-06 5.9e-06 1.4e-05 2.0e-05 2.5e-05 2.9e-05 3.0e-05 3.0e-05];
pos=pos(101:300);
faNorm=faNorm(101:300);

pos=pos*thk/2;
fa=fa*faNorm;


if z>max(pos) || z<min(pos)
    sart = 0;
    return;
end

z=z-min(pos);
pos=pos-min(pos);  % start from 0


if v==0
    
    fa=fa*pi/180;
    fa0=interp1(pos,fa,z);
    
    e1=exp(-TR./T1);
    sart=sin(fa0).*(1-e1)./(1-cos(fa0)*e1);  %% already the transverse magnetization after flip
else
    
    step=v*TR;
    
    dist= z/cos(tilt);
    nRF = ceil(dist./step);
    
    
    da =  dist - (0:nRF-1)*v*TR;% total distance travelled after enterring the slice
    
    za=da*cos(tilt);
    
    
    faa=interp1(pos,fa,za);
    
    faa=faa(end:-1:1);  % flip the order of faa
    
    sz=sz_init;
    
    
    for i=1:nRF-1
        
        sz=sz*cos(faa(i)/180*pi);
        
        sz=1+(sz-1)*exp(-TR/T1);
        
    end
    
    sart = sz*sin(faa(nRF)/180*pi);
    
end

