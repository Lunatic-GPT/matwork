function img=draw_image_roicirc(d,range,show_label,rois)
% rois: an cell array of roi's
if ~exist('show_label','var')
    show_label=true;
end

 id_all=[];
 pos_all=[];
   
 count=0;
d= ri_d1(d);
roi_circles=cell(1,length(rois));
for j=1:length(rois)
    
    roi=ri(rois{j});
    
    if max(roi(:))==1
        roi=clusterize2(roi>0);
    end
    
     
    id=unique(roi(:));
    id(id==0)=[];
    nroi=length(id);
    id_all=[id_all;id];
    roi_circles{j}=0*roi;
   
    pos=zeros(2,nroi);
    
    for i=1:nroi
        count=count+1;
        com=ind2subb(size(roi),find(roi==id(i)));
        pos(:,i)=mean(com,1);
        tmp=mask_circle(size(roi),5,pos(:,i),1,true);
        roi_circles{j}(tmp>0)=1;
    end
    pos_all=cat(2,pos_all,pos);
end

if nargout==0
    draw_image_roi(d,range,roi_circles);
else
    img=draw_image_roi(d,range,roi_circles);
end
if show_label
    for i=1:nroi
        text(pos_all(2,i)-5,pos_all(1,i)-9,num2str(id_all(i)),'FontSize',10,'Color',0*[1,1,1]);
    end
end



