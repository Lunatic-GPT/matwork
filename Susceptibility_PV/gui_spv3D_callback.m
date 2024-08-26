function gui_spv3D_callback(params,sbutton,fit_vessel_phi_switch,rad0,dchi0)

if ~exist('fit_vessel_phi_switch','var')
    fit_vessel_phi_switch=false;
end

if isa(params,'char')  %for bsub
    load(params);
end

init=get(params,'Init. [r,chi]');
init=str2num(init);
if ~exist('rad0','var')
    if isempty(init)
       rad0=0.4;
    else
       rad0=init(1); 
    end
end

if ~exist('dchi0','var')     
    if isempty(init)
       dchi0=0.2;
    else
       dchi0=init(2); 
    end
end


if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end
if strcmp(sbutton,'upload + run')
    sbutton='upload data upload script run';
end

global killdevil;  %matlab version too low on killdevil for line2Distance; 
    
  global skip_fit;
  skip_fit=false;
if ~isempty(strfind(sbutton,'determine_mask_path'))
    
    maskfile=get(params,'mask file');
    
    prefix=strtok(maskfile,'.');
    
    if ~exist(['ThPth_',prefix,'_Curv_Dv.mat'],'file')
        ThinningPathFind(maskfile,false);
        load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
        
        pathCurved_nofix(['ThPth_',maskfile]);  %fixedpath
        fixGap(['ThPth_',prefix,'_Curv']);   %probv
    else
        load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
    end
    %% check path length
    for i=1:length(i_ind_path)
            len(i)=pathLength(i_ind_path{i},ind{i},voxsize,size(c));
            fprintf('Path Length %d: %f\n',i,len(i));
    end
        
    voxsize=ri( ['ThPth_',prefix,'_Curv_Dv.mat'],'','','voxsize');
    
    rad = get(params,'vessel ROI radius (mm)');
    [m,theta0,phi0]=get_fitmask4vroi(['ThPth_',prefix,'_Curv_Dv.mat'],str2num(rad)/mean(voxsize));
    
    save(sprintf('%s_rad%s.mat',prefix,rad),'m','theta0','phi0','len');
    
    disp('determine mask done!');
end

if  ~isempty(strfind(sbutton,'save params'))
    
    resname=par_fname(params,rad0,dchi0);
    
    save(resname,'params');
    
end

if  ~isempty(strfind(sbutton,'upload data'))
    
    [sidd,swi]=fileparts(pwd);
    [pvs,sid]=fileparts(sidd);
    iscan=str2num(swi(5:8));
    
   
   mkdir_remote(sid,killdevil);
   mkdir_remote([sid,'/',swi],killdevil);
   mkdir_remote([sid,'/',sprintf('PHA_IMAGES_%04d',iscan-2)],killdevil);
   mkdir_remote([sid,'/',sprintf('MAG_IMAGES_%04d',iscan-3)],killdevil);

   
   winscp_put(sprintf('../PHA_IMAGES_%04d/TE1.mat',iscan-2),[sid,'/',sprintf('PHA_IMAGES_%04d',iscan-2)],killdevil);
   winscp_put(sprintf('../PHA_IMAGES_%04d/TE2.mat',iscan-2),[sid,'/',sprintf('PHA_IMAGES_%04d',iscan-2)],killdevil);
   winscp_put(sprintf('../MAG_IMAGES_%04d/TE1.mat',iscan-3),[sid,'/',sprintf('MAG_IMAGES_%04d',iscan-3)],killdevil);
   winscp_put(sprintf('../MAG_IMAGES_%04d/TE2.mat',iscan-3),[sid,'/',sprintf('MAG_IMAGES_%04d',iscan-3)],killdevil);
   mask=get(params,'manual mask file');
   winscp_put(mask,[sid,'/',swi],killdevil);
   
end


if  ~isempty(strfind(sbutton,'upload script'))
%     mask=get(params,'manual mask file');

  % winscp_put_mat('gui_spv3D_callback.m');
   
  
    [sidd,swi]=fileparts(pwd);
    [pvs,sid]=fileparts(sidd);
    
   
    gui_spv3D_callback(params,'save params');
    
    
     resname=par_fname(params,rad0,dchi0);
   
    winscp_put(resname,[sid,'/',swi],killdevil);

    
    cmd=sprintf('gui_spv3D_callback(''%s'',''calc_from_pattern'')',resname);
    shname=script_name(params,rad0,dchi0);
  make_matlab_shellscript(killdevil,shname,cmd,8,24);
  winscp_put(shname,[sid,'/',swi],killdevil);
end


if  ~isempty(strfind(sbutton,'run'))
    
    [sidd,swi]=fileparts(pwd); 
    [pvs,sid]=fileparts(sidd);
  
    shname=script_name(params,rad0,dchi0);
    run_cmd_remote(sprintf('source %s',shname),[sid,'/',swi],killdevil);
end
 
if  ~isempty(strfind(sbutton,'get results'))
 [sidd,swi]=fileparts(pwd); 
    [pvs,sid]=fileparts(sidd);
 
    resname=save_name(params,rad0,dchi0);
  
  winscp_get(killdevil,[sid,'/',swi],resname,2);
  
end
if  ~isempty(strfind(sbutton,'calc_from_pattern'))
    
global fit_vessel_phi;
fit_vessel_phi=fit_vessel_phi_switch;
    warning('off','all');
    te4s=get(params,'TE for calc susc');
    
    te4s=str2num(te4s);

    B0=get(params,'B0 (T)');

    for iTE=1:length(te4s)
        
        loc_d=pwd;
        scan=loc_d(end-7:end-4);
        scan=str2num(scan);
        
        fpha=fullfile('..',sprintf('PHA_IMAGES_%04d',scan-2),sprintf('TE%d',te4s(iTE)));
        fmag=fullfile('..',sprintf('MAG_IMAGES_%04d',scan-3),sprintf('TE%d',te4s(iTE)));
        
        
        mag(:,:,:,iTE)=ri_mat_dcm(fmag);
        [ph(:,:,:,iTE),vox_size]=ri_mat_dcm(fpha);
        
     %   TE(iTE)=readdPar(fmag,'EchoTime');
        
    end
    mag=mag/max(mag(:))*100;
    TE =[7.5900,15];
    warning('on','all');
   % mask_file=get(params,'mask file');
    manual_mask_file=get(params,'manual mask file');
    rad=get(params,'vessel ROI radius (mm)');
    %mask_file = sprintf('%s_rad%s.mat',prefix,rad);
    
    mask=ri_d1(manual_mask_file,'','','m');
    theta0=ri(manual_mask_file,'','','theta0');
    phi0=ri(manual_mask_file,'','','phi0');
    len=ri(manual_mask_file,'','','len');
    phrange=max(ph(:))-min(ph(:));
    cmplxd=double(mag).*exp(1i*double(ph)/double(phrange)*2*pi-1i*pi);
        
    lnew = 32;
    
    vessels=get(params,'vessels');
   mask_wm=get(params,'wm mask');
   mask_wm=ri_d1(mask_wm);
   
  [vessels,seg]=strtok(vessels,'-');
 
  vessels=str2num(vessels);
  if ~isempty(seg)
      seg=str2num(seg(2:end));
  end
  
    if isempty(vessels)
        vessels=1:length(theta0);      
    end
      nves=length(vessels);
    cd_fit=zeros(lnew,lnew,lnew,nves,'single');
    cd_data=zeros(lnew,lnew,lnew,nves,'single');
    roi=zeros(lnew,lnew,lnew,nves,'uint16');
    fitres=[];
    nseg=[];
    for i=1:length(vessels)  % pattern search for i>nves
        tic;
        disp(i);
        iv=vessels(i);
        maskv=mask==iv;
        mask2=bwmorph3d(maskv,'thicken',2);
        
        St0=mean_roi(abs(cmplxd),mask2&maskv==0)';
        
       
        
        T2s_v = dchi2T2(dchi0*1e-6)*1000;
        T2s_t = 27; % CSF
        
        par.bgPhase=0;
     
        par.St=St0;
        par.Swm_Ref=mean_roi(abs(cmplxd),mask_wm)';
        Sv0=par.Swm_Ref.*exp(-TE/T2s_v).*exp(TE/T2s_t)*1.05;% 1.05 is the blood-WM partition coefficient
        par.Sv=Sv0(1);
        par.theta=theta0(iv);
        par.phi=phi0(iv);
        par.T2s_v=T2s_v;
        par.T2s_t=T2s_t;
        interp=get(params,'interp');
        par.roi_rad=str2num(rad);
        par.len=len(iv);
        par.fit_S =get(params,'fit S');
        if isempty(seg)
           nseg=floor(par.len/vox_size(1));
           par.seg =  1:nseg;
        else
            par.seg = seg;
        end
        par.cm=[0,0];
      %  nseg=floor(par.len/vox_size(1));

              par.dchi=dchi0;
              par.rad=rad0;
              tmp = calc_vessel_susceptibility_FromPattern3D(cmplxd,maskv,TE,B0,vox_size,vox_size,lnew,interp,par);
            % tmp=0;
              fitres=cat(2,fitres,tmp);
              nseg(i)=length(tmp);
              toc;
              fprintf('Remaining time = %f\n',toc*(nves-i));
    
    end
    
    resname=save_name(params,rad0,dchi0);

    save(resname,'vox_size','TE','fitres','par','nseg');
    
end

% assumed a single initial value condition
if  ~isempty(strfind(sbutton,'show pattern'))
    
    
global fit_vessel_phi;
fit_vessel_phi=fit_vessel_phi_switch;
    warning('off','all');
    te4s=get(params,'TE for calc susc');
    
    te4s=str2num(te4s);

    B0=get(params,'B0 (T)');

    for iTE=1:length(te4s)
        
        loc_d=pwd;
        scan=loc_d(end-7:end-4);
        scan=str2num(scan);
        
        fpha=fullfile('..',sprintf('PHA_IMAGES_%04d',scan-2),sprintf('TE%d',te4s(iTE)));
        fmag=fullfile('..',sprintf('MAG_IMAGES_%04d',scan-3),sprintf('TE%d',te4s(iTE)));
        
        
        mag(:,:,:,iTE)=ri_mat_dcm(fmag);
        [ph(:,:,:,iTE),vox_size]=ri_mat_dcm(fpha);
        
     %   TE(iTE)=readdPar(fmag,'EchoTime');
        
    end
    mag=mag/max(mag(:))*100;
    
    TE =[7.5900,15];
    warning('on','all');
   % mask_file=get(params,'mask file');
    manual_mask_file=get(params,'manual mask file');
    
   % if isempty(manual_mask_file)
    %    manual_mask_file=mask_file;
    %end
    
    %prefix=strtok(manumask_file,'.');
    rad=get(params,'vessel ROI radius (mm)');
    %mask_file = sprintf('%s_rad%s.mat',prefix,rad);
    
    mask=ri_d1(manual_mask_file,'','','m');
    theta0=ri(manual_mask_file,'','','theta0');
    phi0=ri(manual_mask_file,'','','phi0');
  len=ri(manual_mask_file,'','','len');
    phrange=max(ph(:))-min(ph(:));
    cmplxd=double(mag).*exp(1i*double(ph)/double(phrange)*2*pi-1i*pi);
        
    lnew = 32;
    
    vessels=get(params,'vessels');
   
  [vessels,seg]=strtok(vessels,'-');
 
  vessels=str2num(vessels);

      seg=str2num(seg(2:end));
  
   iv =vessels;
        maskv=mask==vessels;
        mask2=bwmorph3d(maskv,'thicken',2);
        
        St0=mean_roi(abs(cmplxd),mask2&maskv==0)';
 
        T2s_v = 1000/275;
        T2s_t = 300;
        Sv0=St0.*exp(-TE/T2s_v).*exp(TE/T2s_t);
       
        par.bgPhase=0;
        par.Sv=Sv0(1);
        par.St=St0;
        par.cm=[0,0];
    debug=1;
    if debug
       tmp=load('Results_roiR_1.8mm_v1-11_r0_0.4_chi0_0.2_interp10_fit_S.mat','fitres');
       par.cm=tmp.fitres.res(3:4);
       par.Sv=tmp.fitres.res(5);
       par.St=tmp.fitres.res(6:7);
       
       
        
    end
        
        par.theta=theta0(iv);
        par.phi=phi0(iv);
         par.T2s_v=T2s_v;
        interp=get(params,'interp');
        par.roi_rad=str2num(rad);
        par.len=len(vessels);
        par.dchi=dchi0;
              par.rad=rad0;

            par.seg = seg;
    par.fit_S = true;
        skip_fit=true;
              tmp = calc_vessel_susceptibility_FromPattern3D(cmplxd,maskv,TE,B0,vox_size,vox_size,lnew,interp,par);
      hfig=get(params,'1st fig#');      
     figure(hfig);
     compare_cd(tmp.cd_rot(:,:,1),tmp.cd_fit_rot(:,:,1),tmp.roi_circ,[-4,13],[-180,180]);
    
      figure(hfig+1);
      compare_cd(tmp.cd_rot(:,:,2),tmp.cd_fit_rot(:,:,2),tmp.roi_circ,[-4,13],[-180,180]);
   
      resid=vec(tmp.cd_rot-tmp.cd_fit_rot,tmp.roi_circ);
      
      resid=[real(resid);imag(resid)];
      resid=sos(resid,1)/sqrt(length(resid));
   
      fprintf('Resid = %f\n',resid);
   
end
if  ~isempty(strfind(sbutton,'show fit'))  %need to modify
   roi_rad=get(params,'vessel ROI radius (mm)');
   vsl=get(params,'vessels');
   
   resname=save_name(params,rad0,dchi0);
   load(resname);

   res=structarray(fitres,'res');
   resid=structarray(fitres,'resid');
   exitflag=structarray(fitres,'exitflag');
   resid=reshape(resid,[length(resid)/length(fitres),length(fitres)]);
   resid=sos(resid,1)/sqrt(length(resid));
   
   h=get(params,'1st fig#');
   
   
   figure(h);
   i=get(params,'seg');
   compare_cd(fitres(i).cd_rot(:,:,1),fitres(i).cd_fit_rot(:,:,1),fitres(i).roi_circ,[-4,13],[-180,180]);
   figure(h+1);
   compare_cd(fitres(i).cd_rot(:,:,2),fitres(i).cd_fit_rot(:,:,2),fitres(i).roi_circ,[-4,13],[-180,180]);
   
   fprintf('Fitting results: ');
   disp(res(i,:));
   fprintf('Residual = %f\n',resid(i));
   fprintf('ExitFlag = %d\n',exitflag(i));
   
   if length(fitres)>1
   figure(h+2);
   subplot(3,1,1);myPlot(1:length(fitres),res(:,1),'ko-','seg','Radius (mm)','');
   subplot(3,1,2);myPlot(1:length(fitres),res(:,2),'ko-','seg','Chi (ppm)','');
   subplot(3,1,3);myPlot(1:length(fitres),resid,'ko-','seg','Residual','');
   end
end


disp('gui_spv3D_callback done');

function [m,theta,phi]=get_fitmask4vroi(Path_mFile,rad)

i_ind_path=ri(Path_mFile,'','','i_ind_path');
c=ri(Path_mFile,'','','c');
sz=size(c);
ind=ri(Path_mFile,'','','ind');

m=zeros(sz,'uint16');
for i=1:length(i_ind_path)
    
    iind=i_ind_path{i}(:);
    sub=ind2subb(sz,ind{i}(iind));
    center=mean(sub,1);
    
    [pca_res,score,latent]=pca(sub(:,:));
    
    len=sos(sub(1,:)-sub(end,:),2);
    
    [theta(i),phi(i)]=unitVec2thetaPhi(pca_res(:,1));
    
    if rad>0
    m_tmp=get_fitmask(theta(i),phi(i),center,len,rad,sz);
    
    m(m_tmp>0)=i;
    else
      m(ind{i}(iind))=i;  
    end
end



function m=get_fitmask(theta,phi,center,len,rad,sz)
% center, len, and rad in units of voxsize; voxel positions 1:n
% radius of the cylindrical mask

[x,y,z]=  meshgrid2(1:sz(1),1:sz(2),1:sz(3));

unitVec=thetaPhi2unitVec(theta,phi);

lin(1,:)=center-unitVec/2*len;
lin(2,:)=center+unitVec/2*len;


[dist,pos]=distance2Line(lin,cat(4,x,y,z));

m=zeros(sz,'uint16');


m(dist<=rad & pos>=0 & pos<=len)=1;



function resname=save_name(params,rad0,dchi0)

    vessels=vessel_name(params);
    interp=get(params,'interp');
    rad=get(params,'vessel ROI radius (mm)');
    resname=sprintf('Results_roiR_%smm%s_r0_%s_chi0_%s_interp%d',num2str(rad),vessels,...
         num2str(rad0),num2str(dchi0),interp);
    global fit_vessel_phi;
    if fit_vessel_phi
     resname=[resname,'_fitphi'];
    end
    
    
    fit_S=get(params,'fit S');
    if fit_S
      resname=[resname,'_fit_S_calcT2v.mat'];    
    else
        resname=[resname,'_calcT2v_wmRef.mat'];    
        
    end
    
    
    

function resname=vessel_name(params)

  vessels=get(params,'vessels');


  vessels=strrep(vessels,' ','_');
    resname=sprintf('_v%s',vessels);




    function shname=script_name(params,rad0,dchi0)
        
        vessels=vessel_name(params);
        fit_S=get(params,'fit S');
        if ~fit_S
            shname=sprintf('spv3D%s_rad0_%s_dchi0_%s.sh',vessels,num2str(rad0),num2str(dchi0));
        else
            shname=sprintf('spv3D%s_rad0_%s_dchi0_%s_fit_S.sh',vessels,num2str(rad0),num2str(dchi0));
        end
        


   function resname=par_fname(params,rad0,dchi0)
             vessels=vessel_name(params);
    fit_S=get(params,'fit S');
    
    interp=get(params,'interp');
    if fit_S
        resname=sprintf('par%s_rad0_%s_chi0_%s_interp%d_fit_S.mat',vessels,num2str(rad0),num2str(dchi0),interp);
    else
        resname=sprintf('par%s_rad0_%s_chi0_%s_interp%d.mat',vessels,num2str(rad0),num2str(dchi0),interp);
    end
    
    







