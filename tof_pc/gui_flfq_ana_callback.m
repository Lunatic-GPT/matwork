function gui_flfq_ana_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

if strcmp(sbutton,'gen vessel mask')
    
    mask_factor=get(params,'vessel mask threshold factor');
    fname=get(params,'phase file');
    wm_mask=get(params,'WM mask');
    
    prefix=fname(1:end-4);
    d=ri(fname);
    m=ri(wm_mask,'','','d');
    csize=get(params,'cluster size');
    roi=mask_threshold(d,m,mask_factor,csize,1,1);
    save([prefix,'_vmask.mat'],'roi');
    
elseif strcmp(sbutton,'clear')
    
    ph_file=get(params,'phase file');
    
    
    dir_name=fileparts(ph_file);
    use_peak=get(params,'max pc pixel');
     if use_peak
        save_name=fullfile(dir_name,'results_retro_time_course_usepeak.mat');
    else
        save_name=fullfile(dir_name,'results_retro_time_course.mat');
     end
    
    if exist(save_name,'file')
        delete(save_name);
    end
        
        
elseif strcmp(sbutton,'Pulsatility time course')
   
        wm_mask=get_fpattern(params,'WM mask');
       vessel_mask=get_fpattern(params,'vessel mask');
       ph_file=get_fpattern(params,'phase file');
    meas=get_fpattern(params,'protocol');
    
    VENC = readsPar(meas,'nVelocity');
    bipolar= readsPar(meas,'alFree[21]');
    if ~isempty(bipolar)
        VENC=VENC/2;
    end
    
    use_peak=get(params,'max pc pixel');
    
    dReadoutFOV=readsPar(meas,'asSlice[0].dReadoutFOV');
    dPhaseFOV=readsPar(meas,'asSlice[0].dPhaseFOV');
    
    lRO=readsPar(meas,'lBaseResolution');
    
    lPE = readsPar(meas,'lPhaseEncodingLines');
    
    vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];
    
    
    %       vox_size=get(params,'vox size (mm)');
    interp_factor=get(params,'interp factor');
    bg_size=get(params, 'bg size');
    
    bg_size_i=round(mean(bg_size./vox_size.*interp_factor));
    num2deg=get(params, 'num2deg');
    
    dir_name=fileparts(ph_file);
    
    if use_peak
        save_name=fullfile(dir_name,'results_retro_time_course_usepeak.mat');
    else
        save_name=fullfile(dir_name,'results_retro_time_course.mat');
    end
    if ~exist(save_name,'file')
        [vall,vall_samebg,vall_nobg,vall_bg]=retro_time_course(ph_file,vessel_mask,wm_mask,VENC,use_peak,ceil(bg_size_i),num2deg);
        %%
        
        dph=ri(ph_file);
        try
            dwm=ri(wm_mask,'','','d');
        catch
            dwm=ri(wm_mask);
        end
        eph=mean_roi(std(double(dph),[],4),dwm);
        ev=eph*VENC/180*num2deg;
        save(save_name, 'vall', 'vall_samebg','vall_nobg','vall_bg', 'ev');
    else
        load(save_name);
    end
   
    v_exclude=get(params,'excluded vessels');
   % v_exclude=0;
   % v_exclude=[1,3,4,8];
    v_exclude(v_exclude==0)=[];
    vall(v_exclude,:)=[];
    vall_samebg(v_exclude,:)=[];
    vall_nobg(v_exclude,:)=[];
    
    %v_exclude=get(params,'excluded vessels');
  
    ind_min=get(params,'ind min');
    ind_max=get(params,'ind peak');
    
    
    %%
    
   % nshift=ceil(size(vall,2)/2)-1;   % the pulse ox peak occurs around the center of the time course
   nshift=get(params,'nshift');
   
    %    nshift=get(params,'time point shifts');
    h2=figure(104);
    [puls_all,vall]= plot_retro_time_course(vall,vall_samebg,vall_nobg,[],h2,nshift,ev,ind_min,ind_max);
   % [puls_all,vall]= plot_retro_time_course(vall_nobg,vall_nobg,vall_nobg,[],h2,nshift,ev,ind_min,ind_max);
    
  %  fprintf('pulsatility = %f\n',puls);
    
    %saveas(h2,fullfile(dir_name,'v_vs_cPhase'));
    savetiff(fullfile(dir_name,'v_vs_cPhase'),h2);
elseif strcmp(sbutton,'flow')
    
    
    mid=get(params,'MID');
    if isempty(mid)     
        Flow_PartialVolume(params);
    else
        interp=get(params,'interp factor');
        for i=1:length(mid)
            dir_str=dir(sprintf('meas_MID%d_*',mid(i)));
            cd(fullfile(dir_str.name,['interp',num2str(interp(1)),'_',num2str(interp(2))]));
            Flow_PartialVolume(params);
            cd('../..');          
        end
        
    end
elseif strcmp(sbutton,'ICA flow')
    
    Flow_PartialVolume_ICA(params);
    
elseif strcmp(sbutton,'show fit')
    
     fname=get_fpattern(params,'results');
     data=ri(fname,'','','data');
     cd_fit=ri(fname,'','','cd_fit');
     for i=1:size(cd_fit,3)
      im_plot=cat(3,real(data(:,:,i)),imag(data(:,:,i)),real(cd_fit(:,:,i)),imag(cd_fit(:,:,i)));
      im_plot=repmat2(im_plot,5);
      figure;imshow4(im_plot,[],[2,2]);
      title(num2str(i));
   
     end
     
     
elseif strcmp(sbutton,'pulsatility (partial volume correction)')
  
    mid=get(params,'MID');
    if isempty(mid)     
        Flow_PartialVolume4Pulsatility(params);
    else
        interp=get(params,'interp factor');
        for i=1:length(mid)
            dir_str=dir(sprintf('meas_MID%d_*',mid(i)));
            cd(fullfile(dir_str.name,['interp',num2str(interp(1)),'_',num2str(interp(2))]));
            Flow_PartialVolume4Pulsatility(params);
            cd('../..');          
        end
        
    end
    
elseif strcmp(sbutton,'save params')
    
    disp('');  % do nothing, always save

elseif strcmp(sbutton,'clear file path')

   set(params,'phase file','');
   set(params,'WM mask','');
   set(params,'vessel mask','');
   set(params,'mean mag file','');
   set(params,'mean phase file','');
   set(params,'data dir','');
   

end

save gui_flfq_ana_params params
disp('gui_flfq_ana_callback done');
