function [contrast,s]=contrast_mask2bg(mask,img,bg,roi)
% mask: the mask 
% img: the image file name
% bg: an mask for the backgroun regions outside the head
% roi: only consider the mask voxels that belong to the roi, if exist;
% [contrast,s]=contrast_mask2bg(mask,img,bg,roi)
% s: a struct containing spvs, sbg,snois,sdnois

T2=ri(img);
roi_bg=ri(bg);    
mask=ri(mask);
roi=ri(roi);
T2=double(T2);
  
ext=m_ext_calc(mask);
  
if exist('roi','var')  %added 1/17/2021
    ext=ext&roi>0;
    mask=mask>0&roi>0;
end    
     
spvs=mean(T2(mask));

sbg=mean(T2(ext>0));
sd_nois=std(T2(roi_bg>0));
snois=mean(T2(roi_bg>0));

s.spvs=spvs;
s.sbg=sbg;
s.snois=snois;
s.sd_nois=sd_nois;
contrast=(spvs-sbg)./snois;
 

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
