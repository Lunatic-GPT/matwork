function calc_Cnorm_air_vs_carb


for i=1:46
  for j=1:2
      disp([i,j]);
      if j==1
          T2=ri(sprintf('../../T2nii/T2_PVS%02d.nii.gz',i));
          pvs=ri(sprintf('../mask_koji_ver2/pvsMask_PVS%02d.nii.gz',i));
          roi=ri(sprintf('../ROIResults_air/ROIs_upsampled/ROIs_PVS%02d_us.nii.gz',i));
          roi_wm=ri(fprintf('../ROIResults_Air/WMROIs/WM_PVS%02d_Carbogen.nii.gz',i));
      else
          
          if ~exist(sprintf('../../T2nii/T2_PVS%02d_CARBOGEN.nii.gz',i),'file')
              continue;
          end
          
          T2=ri(sprintf('../../T2nii/T2_PVS%02d_CARBOGEN.nii.gz',i));
          pvs=ri(sprintf('../mask_koji_ver2/pvsMask_PVS%02d_Carbogen.nii.gz',i));         
          roi=ri(sprintf('ROIs/ROIs_PVS%02d_Carbogen.nii.gz',i));
          roi_wm=ri(fprintf('WMROIs/WM_PVS%02d_Carbogen.nii.gz',i));
      end
      roi_bg=ri('../../T2nii/background.nii.gz');
      ext=m_ext_calc(pvs);
      
     
      for k=1:8
          if k<8
           roi2=roi==k;
          else
            roi2=roi_wm>0;
          end
          spvs(i,k,j)=mean(T2(pvs>0&roi2));
          sbg(i,k,j)=mean(T2(ext>0&roi2));
          snois(i,k,j)=mean(T2(roi_bg>0));
          
      end
      
  end
end
contrast=(spvs-sbg)./snois;

save Cnorm_air_vs_carb_8ROIs contrast spvs sbg snois
      

function c2=m_ext_calc(c)

ind=find(c>0);
l=2;
ind2=zeros(length(ind)*(2*l+1)^3,1);
count=0;
sz=size(c);
for i=1:length(ind)

    ijk=ind2subb(sz,ind(i));

                if ijk(1)-l<1 || ijk(1)+l>sz(1) || ijk(2)-l<1 || ijk(2)+l>sz(2) ||ijk(3)-l<1 || ijk(3)+l>sz(3)    
                    continue;
                end
                
    for i1=-l:l
        for i2=-l:l
            for i3=-l:l
                
                
                count=count+1;
               ind2(count)=ind(i)+i1+i2*sz(1)+i3*(sz(1)*sz(2));
                
            end
        end
    end
    
end
ind2=setdiff(ind2,ind);

ind2(ind2==0)=[];


c2=c*0;
c2(ind2)=1;
