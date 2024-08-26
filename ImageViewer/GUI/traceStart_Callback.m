function traceStart_Callback(hObject, eventdata, handles)
% hObject    handle to traceStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getappdata(handles.figure1,'udata');
roi=getappdata(handles.figure1,'uroi');


if isempty(roi)
    roi=data*0;
end
c=constants;

rad=get(handles.searchRadius,'String');
rad=str2double(rad);


thr_scaled=get(handles.thresholdSlider,'Value');
thr_range=str2num(get(handles.thresholdRange,'String'));
thr=thr_range(1)+(thr_range(2)-thr_range(1))*thr_scaled;
roiVal=str2num(get(handles.ROIVal,'String'));

croplr=str2num(get(handles.croplr,'String'));
cropud=str2num(get(handles.cropud,'String'));

points=getappdata(handles.figure1,'TracePoints');
path=getappdata(hObject,'TracePath');
mask=getappdata(hObject,'TraceMask');

h = drawpoint();
slice=getappdata(gca,'slice');

sel = get(gcf, 'SelectionType');
if strcmpi(sel, 'alt')
    delete(h);
    return;
end
pos=get(h,'Position');
pos=round(pos);
pos(:,1)=pos(:,1)+croplr(1);
pos(:,2)=pos(:,2)+cropud(1);

thr_scaled=get(handles.thresholdSlider,'Value');
thr_range=str2num(get(handles.thresholdRange,'String'));
thr=thr_range(1)+(thr_range(2)-thr_range(1))*thr_scaled;
hypo=get(handles.Hypointense,'Value');
if hypo
    thr=-thr;
end

setappdata(handles.figure1,'isSaved',false);

if size(points,1)>=1
    twopoints=cat(1,points(end,:),[pos,slice]);
    
    path_tmp=searchpath(twopoints,data,rad,~hypo);
    
    if ~isempty(path_tmp)
        path{end+1}=path_tmp;
        mask_tmp=find_mask_near_path(path_tmp,data,thr,rad);
        mask{end+1}=mask_tmp;
        
        points=cat(1,points,[pos,slice]);    
        
        roi=setv_points(roi,mask_tmp,roiVal);
        
        for i=1:length(path)
          roi=setv_points(roi,path{i},c.vpath);
        end
        
        
        setappdata(handles.figure1,'TracePoints',points);
        setappdata(handles.figure1,'TracePath',path);
        setappdata(handles.figure1,'TraceMask',mask);
        setappdata(handles.figure1,'uroi',roi);
       
        showImages(handles);

    else
        fprintf('path not found in %f s\n',timeout);
    end
else
        points=[pos,slice];
        mask = find_mask_near_path(points,data,thr,rad);
         
        
        roi=setv_points(roi,mask,roiVal); 
         
        roi=setv_points(roi,points,c.vpath); 
         
        
        setappdata(handles.figure1,'TracePoints',points);
        setappdata(handles.figure1,'TracePath',{points});
        setappdata(handles.figure1,'TraceMask',{mask});
        
        setappdata(handles.figure1,'uroi',roi);     
        showImages(handles);
end

delete(h);


  
function res=searchpath(points,data,rad,hypo_hyper)
% returns the points on the path including the second terminal point.


roi=line2roi(points(1,:),points(2,:),size(data),rad);

d=data(roi>0);
pos=ind2subb(size(roi),find(roi(:)>0));
[d,pos]=order_data(d,pos,points); % order data so that the points are the first and second elements in d and pos

res=calcPath(pos,d,hypo_hyper);

function [d,pos]=order_data(d,pos,points)

ind1=find(sum(pos==points(1,:),2)==3);
ind2=find(sum(pos==points(2,:),2)==3);

ind=1:size(pos,1);
ind(ind1)=[];
ind(ind2)=[];
ind=[ind1,ind2,ind];
d=d(ind);
pos=pos(ind,:);


function res=calcPath(pos,data,hypo_hyper)
% hypo_hyper: 0 hypo, 1 hyper; 

[s,w]=find_neighbors2(pos,data,hypo_hyper);

sp=sparse(s(:,1)',s(:,2)',w);
sp(size(sp,1)+1:size(sp,2),:)=0;
sp=tril(sp + sp');

[maxlen,pathmax] = graphshortestpath(sp,1,2,'Directed',false);

res=pos(pathmax,:);



function [p,w] = find_neighbors2(pos,data,hypo_hyper)

p=[];
w=[];
for i=1:size(pos,1)
    pos1=pos(i,:);
    for j=i+1:size(pos,1)
        pos2=pos(j,:);
        if ~any(abs(pos1-pos2)>1)
            p(end+1,:)=[i,j];
            if hypo_hyper
            scale=2/(data(i)+data(j));    
            else
            scale=(data(i)+data(j)); %hypo
            end
            w(end+1)= sqrt(sum(abs(pos1-pos2).^2))*scale;
        end
    end
end










