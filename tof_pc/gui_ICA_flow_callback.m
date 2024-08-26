function gui_ICA_flow_callback(params,sbutton)

if ~exist('sbutton','var')
    sbutton=get(gco,'String');
end

if strcmp(sbutton,'apparent flow')
    
    ph_file=get_fpattern(params,'phase file');
    mag_file=get_fpattern(params,'mag file');
    
    fname=get_fpattern(params,'vessel mask');
    
    
    meas=get_fpattern(params,'protocol');
    
    interp=get(params,'interp factor');
    
    nring = 4;
    VENC = readsPar(meas,'nVelocity');
    
    thr_factor=get(params,'threshold factor');
    
    num2deg=get(params,'num2deg');
    dReadoutFOV=readsPar(meas,'asSlice[0].dReadoutFOV');
    dPhaseFOV=readsPar(meas,'asSlice[0].dPhaseFOV');
    
    lRO=readsPar(meas,'lBaseResolution');
    
    lPE = readsPar(meas,'lPhaseEncodingLines');
    
    vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];
    
    vmask=ri_d1(fname);
    ph=ri(ph_file);
    mag=ri(mag_file);
    ph=double(ph);
    mag=double(mag(:,:,:,1));
    vmask=clusterize2(vmask);
    nroi=max(vmask(:));
    
    roi_all=0*vmask;
    flow=zeros(nroi,nring+1);
    v=zeros(nroi,nring+1);
    area=zeros(nroi,nring+1);
    
    for i=1:nroi
        
        mag_full=mean(mag(vmask==i));
        for j=0:nring
            roi_tmp=vmask==i;
            
            mag_mean=mean(mag(roi_tmp));
            for k=1:j
                roi_tmp=imdilate(roi_tmp,[1,1,1;1,1,1;1,1,1]);
                
            end
             roi_tmp= roi_tmp&mag(:,:,:,1)>mag_mean*thr_factor;
            if j==nring
              roi_all(roi_tmp>0) = roi_tmp(roi_tmp>0);            
            end
                
            
%            flow(i,j+1) = sum(ph(roi_tmp).*mag(roi_tmp)/mag_full)*prod(vox_size)/prod(interp)*VENC*10*num2deg/180;
            flow(i,j+1) = sum(ph(roi_tmp))*prod(vox_size)/prod(interp)*VENC*10*num2deg/180;
            v(i,j+1)= mean(ph(roi_tmp))*VENC*10*num2deg/180;
            area(i,j+1)=prod(vox_size)/prod(interp)*sum(roi_tmp(:));       
            fprintf('vessel %d - %d: flow = %4.2f(cm3/s); mean v = %4.1fmm/s; area = %3.1f mm2; r = %3.1f\n',i,j,flow(i,j+1)/1000,v(i,j+1),area(i,j+1),sqrt(area(i,j+1)/pi));
     
        end
    end
    
    
         [int_adjust,snorm]=calc_intensity_adjust_ICA(params,v(:,1)/10);
    
         
    ph=scale2n(ph,100,[-3000,3000]);
    mag=scale2n(mag,100);
    roi_all=bwmorph(roi_all,'remove');
    [im1,cm]=combine_over_under(ph,roi_all,gray(100),[1,0,0],roi_all);
    im2=combine_over_under(mag(:,:,:,1),roi_all,gray(100),[1,0,0],roi_all);
    
    figure(100);imshow4(cat(3,im1',im2'),cm,[1,2]);
    
    fprintf('Int adjust = \n');
    disp(int_adjust);
    
    save ICA_apparent_flow int_adjust flow v area snorm
    
elseif strcmp(sbutton,'make low res')
    
    frecon=get_fpattern(params,'recon file');
    interp_factor=get(params,'lowres interp factor');
    nonzero = get(params,'nonzero fraction');
    fmask_rm=get_fpattern(params,'mask for signal removal');
    d=ri(frecon,'','','im_sb');
    nois=ri(frecon,'','','dsb1_32');
    mid=strtok_no(strtok(filename(frecon),'.'),'_',2);
    if ~isempty(fmask_rm)
    mask_rm = ri_d1(fmask_rm);
    d=setv_roi(d,mask_rm,0);
    
    end
    
    d=mean(d,6);
    
    
    fd=fft2c(d);
    
    sz=size(fd);
    ld=round(sz(1:2).*nonzero);
    
    interp=round(sz(1:2).*nonzero.*interp_factor)./ld;
    
    fd=crop_xz(fd,ld(1),ld(2));
    nois=[];
    d2=ifft2c(fd);
    
     fvessel=get_fpattern(params,'vessel mask');
     interp2=round(interp,2);
     roi_ica=ri_d1(fvessel);
 
     if isempty(fmask_rm)
         do_interp_dim12(d2,interp,mid,nois,['_',num2str(nonzero)],false);
         fmask_wm= fullfile(sprintf('interp%s_%s_%s',num2str(interp2(1)),num2str(interp2(2)),num2str(nonzero)),'mask_wm.mat');
         fmask_ica=fullfile(sprintf('interp%s_%s_%s',num2str(interp2(1)),num2str(interp2(2)),num2str(nonzero)),'mask_ICA.mat');
     else
         
        
           do_interp_dim12(d2,interp,mid,nois,['_',num2str(nonzero),'_sigRm'],false,roi_ica);
         
         fmask_wm= fullfile(sprintf('interp%s_%s_%s_sigRm',num2str(interp2(1)),num2str(interp2(2)),num2str(nonzero)),'mask_wm.mat');
         fmask_ica=fullfile(sprintf('interp%s_%s_%s_sigRm',num2str(interp2(1)),num2str(interp2(2)),num2str(nonzero)),'mask_ICA.mat');
     end
     
     fwm=get_fpattern(params,'WM mask');
     if ~isempty(fwm)
         mwm_us=interp_mask(fwm,sz,nonzero,interp,interp_factor);
         save(fmask_wm,'mwm_us');
     end
     
     if ~isempty(fvessel)
         mv_us=interp_mask(fvessel,sz,nonzero,interp,interp_factor);
         save(fmask_ica,'mv_us');
     end
else
    

end

save gui_flfq_ana_params params
disp('gui_flfq_ana_callback done');


function mv_us=interp_mask(fvessel,sz,nonzero,interp,interp_factor)

     
     m_vessel=ri_d1(fvessel);
     
     
     step=round(size(m_vessel)./round(sz(1:2).*nonzero.*interp));
     mv_us=m_vessel(1:step(1):end,1:step(2):end);
     
     
     
     sz_out=round(sz(1:2).*nonzero.*interp_factor);
     if size(mv_us,1)>sz_out(1)
         mv_us=mv_us(1:sz_out(1),:);
     elseif size(mv_us,1)<sz_out(1)
         mv_us(end+1:sz_out(1),:)=0;
     end
     
     if size(mv_us,2)>sz_out(2)
         mv_us=mv_us(:,1:sz_out(2));
     elseif size(mv_us,2)<sz_out(2)
         mv_us(:,end+1:sz_out(2))=0;
     end
          
     
     