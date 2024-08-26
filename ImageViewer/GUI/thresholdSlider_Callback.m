
function thresholdSlider_Callback(hObject, eventdata, handles)

% --- Executes on slider movement.

% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


data=getappdata(handles.figure1,'udata');
roi=getappdata(handles.figure1,'uroi');

c=constants;

rad=get(handles.searchRadius,'String');
rad=str2double(rad);

hypo=get(handles.Hypointense,'Value');

thr_scaled=get(handles.thresholdSlider,'Value');
thr_range=str2num(get(handles.thresholdRange,'String'));
thr=thr_range(1)+(thr_range(2)-thr_range(1))*thr_scaled;

if hypo
    thr=-thr;
end

roiVal=str2num(get(handles.ROIVal,'String'));

path=getappdata(handles.figure1,'TracePath');
mask=getappdata(handles.figure1,'TraceMask');
for i=1:length(mask)
    
        roi=setv_points(roi,mask{i},0);
        mask_tmp=find_mask_near_path(path{i},data,thr,rad);
        roi=setv_points(roi,mask_tmp,roiVal);
        mask{i}=mask_tmp;
              
end

for i=1:length(path)
   roi=setv_points(roi,path{i},c.vpath);
end
        
        setappdata(handles.figure1,'TraceMask',mask);
        setappdata(handles.figure1,'uroi',roi);
     setappdata(handles.figure1,'isSaved',false);
       
        showImages(handles);
        

