function roi=draw2roi(handles,sz,slice)

roi=zeros(sz(1:2));
noption=length(get(handles.roiShapes,'String'));
for k=1:noption
    j=1;
    while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(slice,j,k)))
        h=findobj(handles.figure1,'Tag',roi_pos_name(slice,j,k));
        if k<6
        api = iptgetapi(h);
         
        roi_pos = api.getPosition();
        else
           roi_pos=get(h,'Position'); 
           
        end
        delete(h);
        ud=get(handles.cropud,'String');
        ud=str2num(ud);
        lr=get(handles.croplr,'String');
        lr=str2num(lr);
        
        roi_pos(:,1)=roi_pos(:,1)+lr(1);
        roi_pos(:,2)=roi_pos(:,2)+ud(1);
        
        if k==1 || k==2
            tmp=zeros(sz(1:2));
            if size(roi_pos,1)==2
                
                dist=sos(diff(roi_pos,1,1),2);
                x=linspace(roi_pos(1,1),roi_pos(2,1),2*ceil(dist));
                y=linspace(roi_pos(1,2),roi_pos(2,2),2*ceil(dist));
                
                for ix=1:length(x)
                    tmp(round(x(ix)),round(y(ix)))=1;
                end
                
            elseif size(roi_pos,1)==1
                
                tmp(round(roi_pos(1,1)),round(roi_pos(1,2)))=1;
            else
                tmp=roipoly(zeros(sz(1:2)),roi_pos(:,2),roi_pos(:,1));  %x, y
                
            end
            
            
        elseif k==3
            
            [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % x, y
            orig=roi_pos(1:2)+roi_pos(3:4)/2;
            tmp=sqrt((xx-orig(2)).^2+(yy-orig(1)).^2)<roi_pos(3)/2;
        elseif k==4
            
            [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % x, y
            orig=roi_pos(1:2)+roi_pos(3:4)/2;
            tmp=sqrt((xx-orig(2)).^2/roi_pos(4)^2+(yy-orig(1)).^2/roi_pos(3)^2)<1/2;
            
        elseif k==5
            
            posx=[roi_pos(1),roi_pos(1)+roi_pos(3),roi_pos(1)+roi_pos(3),roi_pos(1)];
            posy=[roi_pos(2),roi_pos(2),roi_pos(2)+roi_pos(4),roi_pos(2)+roi_pos(4)];
            tmp=roipoly(zeros(sz(1:2)),posy,posx);
        elseif k==6
          %  tmp=line2roi(
           tmp=zeros(sz(1:2));
          for ix=1:size(roi_pos,1)-1
            tmp_tmp=line2roi(roi_pos(ix,:),roi_pos(ix+1,:),sz(1:2),0);
            tmp(tmp_tmp>0)=1;
          end
          
        end
        
        roi(tmp>0)=tmp(tmp>0);
        
        j=j+1;
    end
    
end
