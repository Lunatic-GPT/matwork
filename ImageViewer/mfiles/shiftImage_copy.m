function res=shiftImage_copy(old,new,image2shift,oname,range)

%shift image2shift in the same manner as the shift of old to new.

if ~exist('range','var')
    range=10;
end
old=ri(old);
new=ri(new);
nii=load_untouch_niigz(image2shift);
%oname=filename_append(image2shift,'_shift');
res=true;
for i=-range:range
   for j=-range:range
      
       tmp=circshift(old,[i,j]);
       sel=~isnan(tmp(:))&~isnan(new(:));
       if ~any(new(sel)~=tmp(sel))
           fprintf('shift found - %d %d\n',i,j);
            
            nii.img=circshift(nii.img,[i,j]);
            save_untouch_niigz(nii,oname);
           return;
       end
       
   end
    
end
res=false;
fprintf('Match did not find\n');