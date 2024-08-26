function ROIFlash_Callback(hObject, eventdata, handles)
% hObject    handle to ROIFlash (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%cValue=get(handles.showROI,'Value');
%set(handles.showROI,'Value',~cValue);


 matchOverlay = get(handles.ROIMatchOverlay,'Value');
if ~matchOverlay
 roi_name='uroi';
else
    roi_name='oroi';
end
sl=  getappdata(gca,'slice');  
roi=getappdata(handles.figure1,roi_name);
roi2(:,:,sl)=expand_roi(roi(:,:,sl));
setappdata(handles.figure1,roi_name,roi2);

showImages(handles);


pause(1);
%set(handles.showROI,'Value',cValue);

setappdata(handles.figure1,roi_name,roi);

showImages(handles);


function roi2=expand_roi(roi)

unq=unique(roi(:));
unq(unq==0)=[];
roi2=0*roi;
for i=1:length(unq)
    tmp=bwmorph(roi==unq(i),'dilate',5);
    roi2(tmp)=unq(i);
    
end

