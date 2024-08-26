function res=mask_vessel_pc_tof(pc,tof,m,sig_pc,sig_tof,nmax,csize,rad,max_clust_dist)
% rad: radius in units of pixels

%method 1: the tof mask is defined based on clustered defined with pc
% method 2: pc and tof masks are defined independently, then check for
% overlap

method = 2;
if method ==1
    m_pc=mask_threshold_1tail(pc,m,sig_pc,nmax,csize);
    
    m_pc=clusterize2(m_pc);
    res=m_pc;
    tof=double(tof);
    for i=1:max(m_pc(:))
        
        m_vessel=m_pc==i;
        center=ind2subb(size(m_pc),find(m_pc==i));
        
        center=mean(center,1);
        
        m_circ=mask_circle(size(m_pc),rad,center,1);
        
        m_bg=m_circ>0 & m_pc==0 &m>0;
        
        
        mag_bg=mean(tof(m_bg));
        sd_bg=std(tof(m_bg));
        
        mag_pk = mean(tof(m_vessel));
        
        if mag_pk<mag_bg+sd_bg*sig_tof
            res(m_vessel)=0;
        end
        
    end
    save mask_vessel_pc m_pc;
else
    
    
    
   for i=1:size(pc,3)  % do slice by slice 
    m_pc=mask_threshold_1tail(pc(:,:,i),m(:,:,i),sig_pc,nmax,csize);
    m_tof=mask_threshold_1tail(tof(:,:,i,1),m(:,:,i),sig_tof,nmax,csize);

       m_pc2(:,:,i)=clusterize2(m_pc);
       m_tof2(:,:,i)=clusterize2(m_tof);
       
    res(:,:,i)=clusters_overlap(m_pc2(:,:,i),m_tof2(:,:,i),max_clust_dist);
    
   end
    %m_pc=clusterize2(m_pc);
    save mask_vessel_pc m_pc2;
    save mask_vessel_tof m_tof2;
    
end




