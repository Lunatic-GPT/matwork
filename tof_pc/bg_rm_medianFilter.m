function [a,b_mean]=bg_rm_medianFilter(dname,roi,winSize)
% fl_fq_pc(dname,roi_name)
% roi_name: can be a mask file/matrix or a dicom folder.  In the later, mask will
% % If not given, all voxels will be considered.


a=ri(dname);
b_mean=0*a;

for k=1:size(a,4)    
    for sl=1:size(a,3)
        a_tmp=a(:,:,sl,k);
        b_mean(:,:,sl,k)=median_filter(a_tmp,roi(:,:,sl),winSize);
    end
end

  a=a-b_mean;
  
  roi=repmat(roi,[1,1,1,size(a,4)]);
  a(roi==0)=0;
  
%%

if isa(dname,'char')

prefix=strtok(dname,'.')
    if strcmp(prefix(end-3:end),'.mat')
        prefix=prefix(1:end-4);
    else
        prefix=strtok(prefix,'.');
    end
    
    save([prefix,'_detrend_mdnFltr.mat'],'a');
end


function res=median_filter(d,roi,l)

res=d.*0;

sz=size(d);


for i=floor(l/2)+1:sz(1)-ceil(l/2)+1
    
    if mod(i,100)==0
        %disp(i);
    end
    for j=floor(l/2)+1:sz(2)-ceil(l/2)+1
        
        indi=i-floor(l/2):i+ceil(l/2)-1;
        indj=j-floor(l/2):j+ceil(l/2)-1;
        
            if roi(i,j)==0
                continue;
            end
            
            dtmp=d(indi,indj);
            
           % roitmp=roi(indi,indj);
            
            roitmp=true(size(dtmp)); % use all voxels; including outside roi
            res(i,j)=median(dtmp(roitmp>0));
    end
end





