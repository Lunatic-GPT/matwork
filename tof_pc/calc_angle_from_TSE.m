function th=calc_angle_from_TSE(fname,voxsize,thk)
% fname: name of pvs masks, different pvs should already have different values
% voxsize: in-plane resolution
% thk: slice thickness
% the result is angle in radian
if isa(fname,'char')
roi=ri(fname);
else
    roi=fname;
end

tmp=sum(sum(roi,1),2);
slice=find(squeeze(tmp));
if length(slice)~=3
    error('Something wrong');
end


nroi=max(roi(:));
th=zeros(1,nroi);
for i=1:nroi
    
    cm=zeros(3,2);
    for j=1:3
        pos=ind2subb(size(roi(:,:,1)),find(roi(:,:,slice(j))==i));
        cm(j,:)=mean(pos,1);
        
    end
    
    th(i)=0.5*(atan(sos((cm(1,:)-cm(2,:)).*voxsize)/thk)+atan(sos((cm(2,:)-cm(3,:)).*voxsize)/thk));
 
    
  
    
end
  %disp(costh);

%%




