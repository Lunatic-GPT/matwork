function PA_Angle2SliceNorm(pa_mask,pvs_mask,max_clust_dist)
% nifti file for PVSDir from PVS_orient
% nifti file for PA_mask
% max_clust_dist: unit pixel

nii_pa=load_untouch_niigz(pa_mask);

if ~exist('max_clust_dist','var')
    max_clust_dist=2/nii_pa.hdr.dime.pixdim(2);% 2mm
end
orient=get_orient_from_nii(pa_mask);
norm=orient.rotmat(:,3);

pvs2pa=reslice_dcm(pa_mask,pvs_mask,[],true);

m_pa=nii_pa.img(:,:,1,1);
iab=clusters_match(m_pa,pvs2pa(:,:,1,1),max_clust_dist,false);


sz=size(nii_pa.img);

out=zeros([sz(1:2),1,4]); % 1: angle (PVS dir from mean of whole PVS); 2: angle (PVS dir from each seg);  3: PA index; 4: PVS index

for i=1:size(iab,1)
    
   roi=pvs2pa(:,:,1,1)==iab(i,2);    
   v_whole=getv_roi(pvs2pa(:,:,1,2:4),roi);
   v_seg=getv_roi(pvs2pa(:,:,1,5:7),roi);
   an_whole=angle_bw_2vec(mean(v_whole,1),norm);
   an_seg=angle_bw_2vec(mean(v_seg,1),norm);
      
   if an_seg>90
       an_seg=180-an_seg;
   end
   
   if an_whole>90
       an_whole=180-an_whole;
   end
   
   out(:,:,1,1)=setv_roi(out(:,:,1,1),m_pa==iab(i,1),an_whole);
   out(:,:,1,2)=setv_roi(out(:,:,1,2),m_pa==iab(i,1),an_seg);
   
   out(:,:,1,3:4)=setv_roi(out(:,:,1,3:4),m_pa==iab(i,1),iab(i,:));
   
end

nii_pa.img=out;
save_name=['Angles_maxd_',num2str(max_clust_dist),'_',strtok(filename(pa_mask),'.')];

save_untouch_niigz(nii_pa,save_name);



