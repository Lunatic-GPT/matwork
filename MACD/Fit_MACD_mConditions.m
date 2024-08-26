function Fit_MACD_mConditions(params,save_name)
% Fit_MACD_2Conditions(params,save_name)
%% params should contain the following fields:
% roi_vessel: vessel mask; n*n*1*2
% roi_wm: white matter mask; n*n*1*2;
% voxSize: acquired voxSize; (cm)
% voxSize_interp: voxSize of the images; (cm)
% thk: slice thickness: (cm)
% mag: magnitude images; 5D
% ph: phase images; (value should be in degrees). 5D; the last dim is the
% nphase
% mag_mean:  4D
% ph_mean:
% FA: flip angle;
% roiRad: radius for the circular ROI encompassing all pixels to be included in model image fitting (cm).
% bgRad: background ROI out- and inner- radii for calculating tissue signal (cm)
% interp_factor: interpolation factor for PSF convolution.
% VENC: (cm/s)
% T1: T1 for blood and WM tissue, respectively. (s)
% T2: T2 for blood and WM tissue, respectively. (s)
% TR: (s)
% TE: (s)
% negPhase: whether the flow produce a negative phase
% flow_pattern: 'Laminar' (default),'BluntedParabolic' or 'Plug'
% plot_fit: plot fit results
% skip_fit: skip fitting
% rad0 and v0: radius and flow velocity in mm and cm/s; initial values for
% fitting; also used for calculating the fitting image when skip_fit is
% true;
% lambda: blood-WM partition coefficient
% heart_rate: heart rate in min-1
% weigth_sameSign

T1=params.T1;  %blood, wm.

if ~isfield(params,'T2')
    T2=[28.5,23.5]*0.001;
else
    T2=params.T2; %  %putaman T2* is 17.6 ms
end

neg_phase=params.negPhase; %false;

thk=params.thk; % cm
TR = params.TR; % in s
TE = params.TE; %

FA=params.FA;
if ~isfield(params,'flow_pattern')
    params.flow_pattern='Laminar';
end
[va,sa,s_static]=prep_interp(TR, T1, T2, TE, thk, FA);

%%

    fit_par=params;
    fit_par.va=va;
    fit_par.sa=sa;
    
roi_vessel=params.roi_vessel;

nt=size(params.ph,5);

nroi=max(roi_vessel(:));
v=zeros(nroi,nt);
rad = zeros(nroi,nt);
ratio_r2v=zeros(1,nroi);

center=zeros(2*nt,nroi);
mbg=zeros(2,nroi,nt);
  
cd_fit=[];
data=[];

if neg_phase
    params.ph=-params.ph;
end
if nroi==0
    return;
end
for i=1:nroi
    %%  
    tic;
    ptmp=params;
    for j=1:nt
      ptmp.ph=params.ph(:,:,1,:,j);
      ptmp.mag=params.mag(:,:,1,:,j);
      ptmp.roi_wm=params.roi_wm(:,:,1,j);
      [data(:,:,i,j),roi,mbg(:,i,j),sdata_debug(i)]=calc_cd(ptmp,roi_vessel(:,:,1,j),i,s_static);
    end
    
    
    fit_par.roi=roi;
    fit_par.data=data(:,:,:,1);
    fit_par.sz=size(roi);
    rs=do_fit(fit_par,params.rad0,params.v0,[0,0]);
    
 
    fit_par.data=data;
    
    [rs,exitflag(i),output(i),cd_fit(:,:,i,:)]=do_fit(fit_par,rs.rad,rs.v,[rs.center(1),0,rs.center(2),0]);
    
    center(:,i) = rs.center*10; % in mm
    v(i,:) = rs.v;  % in cm/s
    rad(i,:)=rs.rad; %% in mm
    ratio_r2v(i)=rs.ratio_r2v;
       
    fprintf('Time remaining %4.3f s (%d/%d)\n',toc*(nroi-i),i,nroi);
end

roi_vessel=params.roi_vessel;
params=rmfield(params,{'mag','ph','roi_vessel'});
save(save_name,'v','rad','ratio_r2v','cd_fit','data','center','params','exitflag','output','roi_vessel','sdata_debug','mbg');

function [rs,exitflag,output,cd_fit]=do_fit(fit_par,rad0,v0,center)
% center in cm; rad0 in mm; v0 in cm/s; 
% center dimensions (1,2*nt)

% rs: center in cm;  rad in mm;  v in cm/s

nt=size(fit_par.data,4);
options=prep_options;

% x(1:2) in cm;  x(3) in mm;  x(4) in 10 cm/s
ratio_r2v=0;
if fit_par.ratio_r2v<0 || nt==1
    fitfunc2 = @(x) MBAC_fitfunc_mCondition(x(1:2),x(3:2+nt),x(3+nt:2+nt*2),fit_par);
    if isfield(fit_par,'skip_fit') && fit_par.skip_fit
        res=[center,rad0*ones(1,nt),v0/10*ones(1,nt)];
    else
        x0=[center,rad0*ones(1,nt),v0/10*ones(1,nt)];
        lb=[-0.05,-0.05,0.01*ones(1,nt),zeros(1,nt)];
        ub=[0.05,0.05,0.25*ones(1,nt),0.5*ones(1,nt)];
        [res,tmp,tmp2,exitflag,output]=lsqnonlin(fitfunc2,x0,lb,ub,options);
    end
    center=res(1:2);
    rad=res(3:nt+2);
    v=res(nt+3:2+nt*2)*10;
else
    fitfunc2 = @(x) MBAC_fitfunc_mCondition(reshape(x(1:2*nt),[nt,2]),x(2*nt+1)+x(2*nt+2)*(x(2*nt+3:3*nt+2)*10-x(2*nt+3)),x(2*nt+3:3*nt+2),fit_par); 
    if isfield(fit_par,'skip_fit') && fit_par.skip_fit
        res=[center,rad0,fit_par.ratio_r2v,v0/10*ones(1,nt)];
    else
        x0=[center,rad0,fit_par.ratio_r2v,v0/10*ones(1,nt)];       
        lb=[-0.05,-0.05,-0.05,-0.05,0.01,0,zeros(1,nt)];
        ub=[0.05,0.05,0.05,0.05,0.35,1,0.5*ones(1,nt)];
        [res,tmp,tmp2,exitflag,output]=lsqnonlin(fitfunc2,x0,lb,ub,options);
    end    
    center=res(1:2*nt);
    rad=res(2*nt+1)+res(2*nt+2)*(res(2*nt+3:3*nt+2)*10-res(2*nt+3));
    v=res(2*nt+3:3*nt+2)*10;
    ratio_r2v=res(2*nt+2);
end

rs.rad=rad; %mm
rs.v=v; %cm/s
rs.center=center; %cm
rs.ratio_r2v=ratio_r2v;

[resid,cd_fit]=fitfunc2(res);

if isfield(fit_par,'plot_fit') && fit_par.plot_fit  %check fit
    figure(122);
    nr=ceil(sqrt(nt));
    for t=1:nt
        subplot(nr,nr,t);
        plot_fit(fit_par.data(:,:,1,t),cd_fit(:,:,1,t),roi);
    end
end

for t=1:nt
    fprintf('d = %s mm; v = %s cm/s; \n',num2str(rad(t)*2),num2str(v(t)));
end

function options=prep_options
if datenum(version('-date'))>=736580
    options = optimoptions('lsqnonlin','MaxIter',Inf,'StepTolerance',1e-6,...
        'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,...
        'OptimalityTolerance',0,'Display','iter',...
        'FiniteDifferenceStepSize',0.01);
else
    options = optimoptions('lsqnonlin','MaxIter',Inf,'TolX',1e-6,...
        'MaxFunEvals',Inf,'TolFun',0,'TolPCG',0,'Display','iter','FinDiffRelStep',0.01);
end

function plot_fit(data,cd_fit,roi)

data(~roi)=0;

cd_fit(~roi)=0;
figure; compare_cd(data,cd_fit);
resid=sos(data(:)-cd_fit(:),1).^2;
disp(resid);

function [cd,roi3,mbg,data_debug]=calc_cd(p,roi_vessel,i,s_static)
% cd: the 4th dimension is now cardiac phase

roi2=roi_vessel==i;
p.bgRad=sort(p.bgRad,'descend');
bgRad_i=round(p.bgRad./p.voxSize_interp);

bg_out=mask_circle(size(roi2),bgRad_i(1),roiCOM(roi2),1);
bg_in=mask_circle(size(roi2),bgRad_i(2),roiCOM(roi2),1);

bg=bg_out&~bg_in;

bgroi=bg&p.roi_wm(:,:)>0&roi_vessel==0;

if sum(bgroi(:))<118 % if partly outside wm, do not pose restriction.
    bgroi=bg;%
end

pos=roiCOM(roi2);
roiRad_i=round(mean(p.roiRad./p.voxSize_interp));  % 1 and 2 should be equal

cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;
cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;

cropr_l=pos(1)-bgRad_i(1)-5:pos(1)+bgRad_i(1)+5;
cropc_l=pos(2)-bgRad_i(1)-5:pos(2)+bgRad_i(1)+5;

cropr_l(cropr_l<1)=1;
cropr_l(cropr_l>size(bg,1))=size(bg,1);

cropc_l(cropc_l<1)=1;
cropc_l(cropc_l>size(bg,2))=size(bg,2);

ph0=p.ph;
mag0=p.mag; 

if size(p.ph,4)>1 % phase images
    ph=tof_bg_rm(squeeze(ph0(:,:,1,2*i-1:2*i,:)),bgroi,[],[],false,true);
    mag=tof_bg_rm(squeeze(mag0(:,:,1,2*i-1:2*i,:)),bgroi,[],[],true,true);
else % phase difference images
    
    sz=size(ph0);
    if length(sz)==2
        sz(3)=1;
    end
    
    nt=size(ph0,5);
    ph=reshape(ph0,[sz(1:3),nt]);   
    ph=tof_bg_rm(ph,bgroi,[],[],false,true);
    mag=zeros([sz(1:2),2,nt]);
    
    magtmp=reshape(mag0(:,:,:,1,:),[sz(1:3),nt]);
    mag(:,:,1,:)=tof_bg_rm(magtmp,bgroi,[],[],true,true);
    
    magtmp=reshape(mag0(:,:,:,2,:),[sz(1:3),nt]);
    mag(:,:,2,:)=tof_bg_rm(magtmp,bgroi,[],[],true,true);
    
end

%mbg=mean_roi(mean(mag,4),bgroi>0);
mbg=mean_roi(mean(mag,3),bgroi>0);
snorm=mbg/s_static*p.lambda;

mag1=mag(:,:,1,:);
mag2=mag(:,:,2,:);

if size(ph,3)==1
    cd=mag2.*exp(1i*ph/180*pi) - mag1;
else
    cd=mag2.*exp(-1i*ph(:,:,2,:)/180*pi) - mag1.*exp(-1i*ph(:,:,1,:)/180*pi);
end

cd=cd/snorm;

%         z(1)=mean(cd(roi2));
%         z(2)=mean(mag1(roi2))/snorm;
%

%% debug purpose
data_debug.ph_crop=ph(cropr_l,cropc_l,:,:,:);
data_debug.mag1_crop=mag(cropr_l,cropc_l,1,:);
data_debug.mag2_crop=mag(cropr_l,cropc_l,2,:);
data_debug.cd_crop=cd(cropr_l,cropc_l,:,:);
data_debug.bgroi_crop=bgroi(cropr_l,cropc_l);
roi3=mask_circle(size(cd),roiRad_i,pos,1);
data_debug.roi_crop=roi3(cropr_l,cropc_l);

roi3=roi3(cropr,cropc);
cd=cd(cropr,cropc,:,:);

