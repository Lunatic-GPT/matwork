function roi=mask_PAinPC(mag,ph,wm_mask,alpha_2tail,voxsize)
% mask_PAinPC_new(mag,phase,wm_mask,alpha_2tail)
% calcu a PA mask based on a series of phase (in radian) and mag images
% tSNR calculated to determine threshold for roi determination.
% par should contain a parameter named alpha_2tail 

if ~exist('alpha_2tail','var')
    alpha_2tail=0.1;
end


par.SNR_mag=mean(mag,4)./std(mag,[],4);

par.patchSize=round(10/voxsize(1));

par.wm_mask=clusterize2(ri_d1(wm_mask),100);

par.data=mean(mag(:,:,:,1:2:end),4);

par.ispc=false;
[mgmn_p22,mgmn_p05]=calc_2masks(par,alpha_2tail);
    
par.data=mean(ph,4);
par.ispc=true;
[phmn_p22,phmn_p05]=calc_2masks(par,alpha_2tail);
  
roi=find_overlap(phmn_p22,mgmn_p22);
      
roi=clusterize2(roi);

save temp roi mgmn_p22 phmn_p22  mgmn_p05  phmn_p05
         
%roi=repmat(roi,[1,1,1,2]); %to display correctly in ITKSnap; no need, can
%be correctly displayed as mask, but not as additional images




function [msqrt,m2]=calc_2masks(par,alpha)
    par.alpha_2tail=sqrt(alpha);
    msqrt=vmask_medianFilter(par);
    par.alpha_2tail=alpha;
    m2=vmask_medianFilter(par);
    m2=clusterize2(m2);
    msqrt=clusterize2(msqrt);
    
function m_ph=find_overlap(m_ph,m_mg)
  m_ph=clusterize2(m_ph);
  m_mg=clusterize2(m_mg);
      
  [m_ph,m_mag]=clusters_overlap(m_ph,m_mg);
    


    
    
    
