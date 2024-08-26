
% --- Executes on button press in ApplyROI.
function ApplyROI_Callback(hObject, eventdata, handles)
% hObject    handle to ApplyROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~exist('dosave','var')
    dosave=true;
end

matchOverlay = get(handles.ROIMatchOverlay,'Value');

if matchOverlay
    roi=getappdata(handles.figure1,'oroi');
else
    roi=getappdata(handles.figure1,'uroi');
end

if matchOverlay
    d=getappdata(handles.figure1,'odata_rs');
    if isempty(d)
        d=getappdata(handles.figure1,'udata');
        set(handles.ROIMatchOverlay,'Value',false);
    end
    
    ind_u2o=getappdata(handles.figure1,'ind_u2o');
else
    d=getappdata(handles.figure1,'udata');
    ind_u2o=0;
end
sz=size(d);
if length(sz)==2
    sz(3)=1;
end
if isempty(ind_u2o)
    ind_u2o=0;
end

val = str2num(get(handles.ROIVal,'String'));

roi_new=saveROI_Callback(hObject, eventdata, handles,false);
applyTo=get(handles.listbox7,'Value');

if applyTo ==1
     
if matchOverlay
    setappdata(handles.figure1,'oroi',roi_new);
else
    setappdata(handles.figure1,'uroi',roi_new);
end

else
    if ~isempty(roi)
      d= setv_roi(d,roi_new>0&roi==0,val);
    else
      d= setv_roi(d,roi_new>0,val); 
    end
    if matchOverlay
        setappdata(handles.figure1,'odata_rs',d);
    else
        setappdata(handles.figure1,'udata',d);
    end
   
end
save ApplyROI_temp_after roi
showImages(handles);