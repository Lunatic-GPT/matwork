 function MouseButtonDown_Callback(hObject, eventdata)
 
 handles=guidata(hObject);
 if strcmp(hObject.SelectionType,'normal') %left button
     
     %handles=getappdata(hObject,'handles');    
     val=get(handles.showOverlay,'Value');
     set(handles.showOverlay,'Value',~boolean(val));
     showImages(handles);
     
 elseif   strcmp(hObject.SelectionType,'extend') %middle button
    
      
 else
     % handles=getappdata(hObject,'handles'); %right button
     val=get(handles.showROI,'Value');
     set(handles.showROI,'Value',~boolean(val));
     showImages(handles);
 end
   
        %     figure(gh);
             