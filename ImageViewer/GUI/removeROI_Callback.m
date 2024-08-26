
% --- Executes on button press in removeROI.
function removeROI_Callback(hObject, eventdata, handles)
% hObject    handle to removeROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%set(handles.figure1,'WindowButtonDownFcn','');
Draw1_Callback(hObject, eventdata, handles);


d=getappdata(handles.figure1,'udata');
sz=size(d);
slice=getappdata(gca,'slice');
m=draw2roi(handles,sz,slice);

ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');




all_slice = false;  % use a gui later;

if ROIMatchOverlay
    roi2=getappdata(handles.figure1,'oroi');
else
    roi2=getappdata(handles.figure1,'uroi');
end

if all_slice
   slice=1:size(roi2,3);
end    
    if ~isempty(roi2)
        roi2(:,:,slice)=double(roi2(:,:,slice)).*repmat(double(m==0),[1,1,length(slice)]);
    end

    
    

if ROIMatchOverlay
    setappdata(handles.figure1,'oroi',roi2);
else
    setappdata(handles.figure1,'uroi',roi2);
end

showImages(handles);