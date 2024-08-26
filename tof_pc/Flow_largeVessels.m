function flow = Flow_largeVessels(meas,sub_dir,mask_vessel,mask_vessel_pvf1,interp_factor,neg_phase)
   

    VENC = readsPar([meas,'.pro'],'nVelocity');
    bipolar= readsPar([meas,'.pro'],'alFree[21]');
    dReadoutFOV=readsPar([meas,'.pro'],'asSlice[0].dReadoutFOV');
    dPhaseFOV=readsPar([meas,'.pro'],'asSlice[0].dPhaseFOV');
    lRO=readsPar([meas,'.pro'],'lBaseResolution');
    lPE = readsPar([meas,'.pro'],'lPhaseEncodingLines');
    vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];
    vox_size=prod(vox_size./interp_factor);   
    if ~isempty(bipolar)
        VENC=VENC/2;
    end
    
    if ~exist('neg_phase','var')
        neg_phase=true;
    end
    cd(meas);
    mid1=strtok_no(meas,'_',2);
    ph_file=sprintf('%s/Phase_%s.mat',sub_dir,mid1);
    mag_file=sprintf('%s/Mag_%s.mat',sub_dir,mid1);
    
    mag=ri(mag_file);
    ph=ri(ph_file);

    if neg_phase
       ph=-ph;
    end
    mag=double(mag);
    ph=double(ph);
    try
    roi=ri([sub_dir,'/',mask_vessel],'','','d');
    catch
    roi=ri([sub_dir,'/',mask_vessel]);
        
    end
    if max(roi(:))==1
        roi=clusterize2(roi);
    end
    
    for i=1:max(roi(:))

    %%
    
    bg=roiCOMNeighbor(roi==i,15*interp_factor);
    
    bgroi=bg&roi==0;
    fprintf('%d bg voxels = %d \n',i,sum(bgroi(:)));
    
    
  %  mbg=bg4medianFilter(mag(:,:,:,1),33*interp_factor,bgroi);
      
   
        
     phbg=bg4medianFilter(ph,33*interp_factor,bgroi);
    
     [mag2,bgmean]=tof_bg_rm_medianFilter(mag(:,:,:,1),roi==i|bgroi,33*interp_factor);
     
     mag2(roi==i|bgroi)=mag2(roi==i|bgroi)+bgmean;
     
    %%
    
    nv_full=2;  % the square size with no partial volume effect
    
    
    %roi2=find_square_no_partial_volume(ph,roi==i,nv_full);
    roi2=ri([sub_dir,'/',mask_vessel_pvf1]);
    
    mag_full=mean(mag2(roi2>0));
    
    flow(i)=sum(mag2(roi==i)./mag_full.*(ph(roi==i)-phbg)*VENC/18000)*vox_size*10;
    
    
    end
    
    
    save flowLargeVessels flow ;
    cd('..');
    
    function roi2=find_square_no_partial_volume(d,roi,nv)
        
        for i=1:size(roi,1)-nv+1
          
            for j=1:size(roi,2)-nv+1
                
               if any(vec(roi(i:i+nv-1,j:j+nv-1))==0) 
                   continue;
               end
                
               if ~exist('i_max','var')
                   i_max=i;
                   j_max=j;
                   d_max=mean(vec(d(i:i+nv-1,j:j+nv-1)));
               else
                  
                   if mean(vec(d(i:i+nv-1,j:j+nv-1)))>d_max
                   
                       i_max=i;
                       j_max=j;
                       d_max=mean(vec(d(i:i+nv-1,j:j+nv-1)));
                   
                   end
               end
                   
                
            end
            
        end
        
        
        
        
        roi2=0*roi;
        
        roi2(i_max:i_max+nv-1,j_max:j_max+nv-1)=1;
        
           
        
        
        
    
    
    
    
    
