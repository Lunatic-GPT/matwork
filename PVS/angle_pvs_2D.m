function [th_pvs,ipvs_all]=angle_pvs_2D(pvs,pvs_dir,slice_norm)    

% th_pvs: angles to the slice normal direction
% ipvs_all: pvs label;

   pvs=ri(pvs);
   
   ipvs_all=unique(pvs(:));
   ipvs_all(ipvs_all==0)=[];
   
    npvs=length(ipvs_all);
    
    
    pvs_dir=ri_d1(pvs_dir);
    pvs_dir(isnan(pvs_dir))=0;
     
    th_pvs=zeros(npvs,1);
    
    for j=1:npvs
        dir_mean=mean_roi(pvs_dir,pvs==ipvs_all(j));
        th_pvs(j)=acos(sum(dir_mean(:).*slice_norm,1))/sos(dir_mean(:),1)*180/pi;    
    end
    