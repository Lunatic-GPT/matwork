function PVS_region_check(pvs, roi,prefix)
% pvs_seg_fname: output from pvs_quantify
% roi: 4d matrix containing size(roi,4) rois or 3d matrix with values of 1:nroi

%pvsind: PVS indices; 1*nroi cell
%pvsind_noamb: PVS indices excluding ambiguous PVS; 1*nroi cell
%roi_label: pvs mask image where the value of each pvs voxel is set to the
%roi index
%roi_nvox: 1*nroi - number of voxel in each roi

if size(roi,4)==1
    roi2=zeros([size(roi),max(roi(:))]);
    
   for i=1:max(roi(:))
     roi2(:,:,:,i)=setv_roi(roi2(:,:,:,i),roi==i,1);
      
   end
    
   roi=roi2;
end


[pvsind,pvsind_noamb,roi_label,roi_nvox] = get_regions(pvs,roi);

save(prefix,'pvsind','pvsind_noamb','roi_label','roi_nvox');




% else
%   fname=filename_append(pvs_seg_fname,suffix);
%     
%  mat_addv(fname,'roi_nvox',roi_nvox);
% end

function res=calc_roi_nvox(roi) 

for i=1:size(roi,4)
   res(i)=sum(vec(roi(:,:,:,i)>0));
end


function  [ind_n,ind_n2,roi_label,roi_nvox] = get_regions(pvs,roi)



roi_nvox=calc_roi_nvox(roi);


nroi=size(roi,4);
ind_n=cell(1,nroi);  % PVS indices 
ind_n2=cell(1,nroi);  % PVS indices excluding ambiguous PVS
roi_label = 0*roi(:,:,:,1);
sz=size(roi);
roi=reshape(roi,[prod(sz(1:3)),sz(4)]);


for i=1:length(pvs)
 
    indp=pvs(i).ind(1:pvs(i).nvox_path);
    
    ntmp=sum(roi(indp,:),1);
    
    if ~any(ntmp>0)
        continue;
    elseif sum(ntmp>0)==1  % PVS only intersect one region
        ind_n{ntmp>0}(end+1)=i;
        ind_n2{ntmp>0}(end+1)=i;
        roi_label(pvs(i).ind)=find(ntmp>0);
    else % PVS intersect more than one region
        
        [nmax,ireg]=max(ntmp);  
        ind_n{ireg}(end+1)=i;
        roi_label(pvs(i).ind)=ireg;
    end  
    
end

disp(' ');


