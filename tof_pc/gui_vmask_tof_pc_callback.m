function gui_vmask_tof_pc_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end
if strcmp(sbutton,'Phase detrend')  || strcmp(sbutton,'Phase+Mag detrend') || strcmp(sbutton,'detrend+int. thr.') || strcmp(sbutton,'detrend+Frangi')
    
    f_pc=get_fpattern(params,'phase file');
    f_mag=get_fpattern(params,'mag file');
    f_mask=get_fpattern(params,'brain mask');
    
    
    %f_pc=dicomread(f_pc);
    %f_mag=dicomread(f_mag);
    
    
    if isempty(f_mask)
        mask_factor=get(params,'brain mask threshold factor');
        fl_fq_retro_pc(f_pc,mask_factor,f_mag);
    else
        try
            f_mask=ri(f_mask,'','','d');
        catch
            f_mask=ri(f_mask);
        end
        
        flipxy=get_default(params,'mask: flip xy',false);
        
        if flipxy
            f_mask=permute(f_mask,[2,1,3,4]);
        end
        
        fl_fq_retro_pc(f_pc,[],f_mask);
    end
end

if strcmp(sbutton,'Mag detrend') || strcmp(sbutton,'Phase+Mag detrend') || strcmp(sbutton,'detrend+int. thr.')|| strcmp(sbutton,'detrend+Frangi')
    f_mag=get_fpattern(params,'mag file');
    
    % prefix=strtok(f_mag,'.');
    
    %    if ~isempty(strfind(f_mag,'.mat'))
    %      prefix=f_mag(1:end-4);
    %    end
    prefix=strtok2(f_mag,'.');
    mag=ri_d1(f_mag,1);
    mag=mag(:,:,:,1);
    
    f_mask=get_fpattern(params,'brain mask');
    
    
    if isempty(f_mask)
        mask_factor=get(params,'brain mask threshold factor');
        
        tof_bg_rm(mag,[],mask_factor,prefix);
    else
        try
            f_mask=ri(f_mask,'','','d');
        catch
            f_mask=ri(f_mask);
        end
        
        flipxy=get_default(params,'mask: flip xy',false);
        
        if flipxy
            f_mask=permute(f_mask,[2,1,3,4]);
        end
        
        tof_bg_rm(mag,f_mask,[],prefix);
    end
    
    
end

if strcmp(sbutton,'Frangi') || strcmp(sbutton,'detrend+Frangi')
    
    f_m=get_fpattern(params,'WM mask');
    try
    mask=ri(f_m,'','','d');
    catch
    mask=ri(f_m);
        
    end
    f_pc=get_fpattern(params,'phase file');
    
    f_mag=get_fpattern(params,'mag file');
    use_detrend=get(params,'Use _detrend');
    if use_detrend
        
        f_pc=strtok2(f_pc,'.');
        f_pc=[f_pc,'_detrend.nii'];
        
        f_mag=strtok2(f_mag,'.');
        f_mag=[f_mag,'_detrend.nii'];
    end
    
    %%
    % .FrangiScaleRange : The range of sigmas used, default [1 8]
    %       .FrangiScaleRatio : Step size between sigmas, default 2
    %       .FrangiBetaOne : Frangi correction constant, default 0.5
    %       .FrangiBetaTwo : Frangi correction constant, default 15
    %       .BlackWhite : Detect black ridges (default) set to true, for
    %                       white ridges set to false.
    %       .verbose : Show debug information, default true
    
    
    interp_factor=get_default(params,'interp factor',1);
    meas=get_fpattern(params,'protocol');
    if ~isempty(meas)
        dReadoutFOV=readsPar(meas,'asSlice[0].dReadoutFOV');
        dPhaseFOV=readsPar(meas,'asSlice[0].dPhaseFOV');
        lRO=readsPar(meas,'lBaseResolution');
        lPE = readsPar(meas,'lPhaseEncodingLines');
        vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];
    else
        vox_size=get(params,'voxel size');
    end
    
    pc=double(ri(f_pc));
    tof=double(ri(f_mag));
    tof=tof(:,:,:,1);
    options.FrangiScaleRange= [0.25,1];
    options.FrangiScaleRatio = 0.25;
    options.BlackWhite=false;
    
    if 0
    v_pc=FrangiFilter2D(pc,options,[],mask);
    thr_ph=get(params,'vesselness threshold (ph)');
    m_pc=v_pc>thr_ph;
    m_pc(mask==0)=0;
    save('vesselness_pc','v_pc');
    
    else
     areamax=get(params,'max vessel area (mm2)');
    
    
    
    csize=get(params,'cluster size thr (mm2)');
    csize_i=round(csize/prod(vox_size)*interp_factor(1)*interp_factor(2));
    nmax=areamax/prod(vox_size)*interp_factor(1)*interp_factor(2);
    nmax=round(nmax);
    
    sig_pc=get(params,'phase threshold factor'); % 95% one -tail
    m_pc=mask_threshold_1tail(pc,mask,sig_pc,nmax,csize_i);
    
    end
    
    dir_str=dir('vesselness_tof.mat');
    
    if isempty(dir_str)
    v_tof=FrangiFilter2D(tof(:,:,:,1),options,[],mask);
    
    save('vesselness_tof','v_tof');
    else
       load('vesselness_tof.mat'); 
    end
    thr_mag=get(params,'vesselness threshold (mag)');
    m_tof=v_tof>thr_mag;
    m_tof(mask==0)=0;
    
    maxsep=get(params,'max overlap separation (mm)');
    maxsep_i=min(round(maxsep./vox_size));
    
    roi=clusters_overlap(m_pc,m_tof,maxsep_i);
    
    
    
    roi=clusterize2(roi);
    
    
    
    rad_roi=get(params, 'roi circle radius (mm)');
    rad_roi_i=round(mean(rad_roi./vox_size.*interp_factor));
    
    roi_circ=0*roi;
    shape_thr=get(params,'shape threshold');
    roi_total=0;
    for i=1:max(roi(:))
        
        mi=roi==i;
        
        dist=cluster_maxdist(mi);
        
        dia=2*sqrt(sum(mi(:))/pi);
        
        if dist/dia>shape_thr
            roi(mi)=0;
            continue;
        end
        
        cm=ind2subb(size(roi),find(mi));
        cm=mean(cm,1);
        roi_circ=roi_circ|mask_circle(size(roi),rad_roi_i,cm,1);
        roi_total=roi_total+1;
    end
    
    %%
    % roi_draw=m;
    % roi_draw(roi_circ>0)=2;
    %draw_image_roi(tof',[0 0.1*max(tof(:))],mask',roi_circ');
    icrop=get(params,'crop(u d l r)');
    
   [img1,cm]= draw_image_roi(v_tof',[0 thr_mag],mask',roi_circ');
   
   img_tof=draw_image_roi(tof',[0,6000],mask',roi_circ');
   
   draw_image_roi_onIMG(tof,[0,6000],v_tof,[0 thr_mag],roi);
   
    figure(101);imshow(cat(2,img1(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4)),img_tof(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4))),cm);
    
    img2=draw_image_roi(pc',[-30,30],mask',roi_circ');
    
    figure(102);imshow(cat(2,img1(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4)),img2(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4))),cm);
    
   [img1,cm]= draw_image_roi(v_tof',[0 thr_mag],mask',mask'& m_tof',mask'& m_pc',mask'& m_tof'&m_pc',roi_circ');
    
    img2=draw_image_roi(pc',[-30,30],mask',mask'&m_tof',mask'& m_pc',mask'& m_tof'&m_pc',roi_circ');
    
    
    
    figure(103);imshow(cat(2,img1(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4)),img2(icrop(1)+1:end-icrop(2),icrop(3):end-icrop(4))),cm);
    
    %%
    fprintf('%d ROIs detected\n',roi_total);
    
    prefix_pc=f_pc(1:end-4);
    %  prefix_mag=f_mag(1:end-4);
    save([prefix_pc,'_Frangi_vmask.mat'],'roi');
    rad_str=strrep(num2str(rad_roi),'.','_');
    save([prefix_pc,'_Frangi_vcircmask_r',rad_str,'.mat'],'roi_circ');
    
    
    
    
    %%
    
end
if strcmp(sbutton,'int. thr.') || strcmp(sbutton,'detrend+int. thr.') 
    f_pc0=get_fpattern(params,'phase file');
    f_pc=f_pc0;
    f_mag=get_fpattern(params,'mag file');
    
    use_detrend=get(params,'Use _detrend');
    if use_detrend
        
        f_pc=strtok2(f_pc,'.');
        f_mag=strtok2(f_mag,'.');
        
        f_pc=[f_pc,'_detrend.mat'];
        f_mag=[f_mag,'_detrend.mat'];
        
    end
    
    sig_pc=get(params,'phase threshold factor'); % 95% one -tail
    sig_tof=get(params,'mag threshold factor'); % 95% one -tail
    
    f_wm=get_default(params,'WM mask',get_fpattern(params,'brain mask'));
    
    
    try
        m=ri(f_wm,'','','d');
    catch
        m=ri(f_wm);
    end
    
    try
        flipxy=get(params,'mask: flip xy');
    catch
        flipxy = false;
    end
    
    if flipxy
        m=permute(m,[2,1,3,4]);
    end
    
    try
        pc=ri(f_pc,'','','d');
    catch
        pc=ri(f_pc);
    end
    tof=ri(f_mag);
    
    csize=get(params,'cluster size thr (mm2)');
    
    meas=get_fpattern(params,'protocol');
    if ~isempty(meas)
        dReadoutFOV=readsPar(meas,'asSlice[0].dReadoutFOV');
        dPhaseFOV=readsPar(meas,'asSlice[0].dPhaseFOV');
        lRO=readsPar(meas,'lBaseResolution');
        lPE = readsPar(meas,'lPhaseEncodingLines');
        vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];
    else
        vox_size=get(params,'voxel size');
    end
    
    
    interp_factor=get_default(params,'interp factor',1);
    if length(interp_factor)==1
        interp_factor=interp_factor*[1,1];
    end
    rad=get(params, 'bg circle radius (mm)');
    
    areamax=get(params,'max vessel area (mm2)');
    
    
    
    rad_i=round(mean(rad./vox_size.*interp_factor));
    csize_i=round(csize/prod(vox_size)*interp_factor(1)*interp_factor(2));
    
    nmax=areamax/prod(vox_size)*interp_factor(1)*interp_factor(2);
    nmax=round(nmax);
    maxsep=get(params,'max overlap separation (mm)');
    
    maxsep_i=min(round(maxsep./vox_size));
    
    roi=mask_vessel_pc_tof(pc,tof,m,sig_pc,sig_tof,nmax,csize_i,rad_i,maxsep_i);
    roi=clusterize2(roi);
    
    rad_roi=get(params, 'roi circle radius (mm)');
    rad_roi_i=round(mean(rad_roi./vox_size.*interp_factor));
    
    roi_circ=0*roi;
    shape_thr=get(params,'shape threshold');
    roi_total=0;
    for i=1:max(roi(:))
        
        mi=roi==i;
        
        dist=cluster_maxdist(mi);
        
        dia=2*sqrt(sum(mi(:))/pi);
        
        if dist/dia>shape_thr
            roi(mi)=0;
            continue;
        end
        
        cm=ind2subb(size(roi),find(mi));
        cm=mean(cm,1);
        roi_circ=roi_circ|mask_circle(size(roi),rad_roi_i,cm,1);
        roi_total=roi_total+1;
    end
    
    %%
    % roi_draw=m;
    % roi_draw(roi_circ>0)=2;
    
    
%    draw_image_roi_onIMG(pc,[-30,30],tof,[0 0.1*max(tof(:))],roi,m);
    
%     draw_image_roi(tof',[0 0.1*max(tof(:))],m',roi_circ');
%     
%     draw_image_roi(pc',[-30,30],m',roi_circ');
    
    
    %%
    fprintf('%d ROIs detected\n',roi_total);
    
    save_data(f_pc0,roi>0,roi_circ>0);
    
    
    
end
% save([prefix_pc,'_vmask.mat'],'m_pc');
% save([prefix_mag,'_vmask.mat'],'m_tof');
if strcmp(sbutton,'calc flow')
    
    
    f_pc=get_fpattern(params,'phase file');
    
    use_detrend=get(params,'Use _detrend');
    if use_detrend
        
        f_pc=strtok2(f_pc,'.');
        
        f_pc=[f_pc,'_detrend.mat'];
        
    end
    
    prefix_pc=strtok2(f_pc,'.');
    
    fname=[prefix_pc,'_tof_vmask.mat'];
    vmask=get_default(params,'vessel mask',fname);
    voxsize=get(params,'voxel size');
    VENC=get(params,'VENC');
    try
        pc=ri(f_pc,'','','d');
    catch
        pc=ri(f_pc);
    end
    vmask=ri(vmask);
    
    vmask=clusterize2(vmask,1);
    nc=max(vmask(:));
    flow=zeros(1,nc);
    for i=1:nc
        flow(i)=sum(pc(vmask==i))*prod(voxsize)*VENC/180*10; %mm3/s
    end
    disp('Flow (mm3/s)');
    disp(flow);
    [dname,pc_name]=fileparts(f_pc);
    save(fullfile(dname,sprintf('flow_%s.mat',pc_name)),'flow');
end
if strcmp(sbutton,'clear file path')
    set(params,'mag file','');
    set(params,'phase file','');
    set(params,'brain mask','');
    try
        set(params,'WM mask','');
    catch
    end
end

    save gui_vmask_tof_pc_params params;

disp([mfilename,' Done']);


function save_data(f_pc,roi,roi_circ)
    
    prefix_pc=f_pc(1:end-4);
    fname_roi=[prefix_pc,'_vmask.mat'];
    fname_roi_circ=[prefix_pc,'_vcircmask.mat'];
    
    if strcmp(f_pc(end-3:end),'.mat')
    
    save(fname_roi,'roi');
    save(fname_roi_circ,'roi_circ');
    
    else
    
        nii=load_untouch_nii(f_pc);
        nii.img=roi;
        save_untouch_nii(nii,fname_roi);
        
        nii.img=roi_circ;
        save_untouch_nii(nii,fname_roi_circ);
        
    end
    
    
    