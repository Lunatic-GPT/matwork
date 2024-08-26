 function WindowKeyPressFcn_Callback(hObject, eventdata)
 if isempty(eventdata.Character)
     return;
 end
 
 handles=guidata(hObject);
 c=constants;
 if eventdata.Character=='a' %add a new point
    traceStart_Callback(hObject, eventdata, handles);   
 elseif eventdata.Character=='z' %draw ROI
   Draw1_Callback(hObject, eventdata, handles);   
 elseif eventdata.Character=='x' %apply ROI
    ApplyROI_Callback(hObject, eventdata, handles);
    
 elseif eventdata.Character=='d' %delete a point and associated path
     points=getappdata(hObject,'TracePoints');
     path=getappdata(hObject,'TracePath'); % path contains a list of points between two seeds, including the second seed. 
                                           % the first element only contain the very first seed.
     mask=getappdata(hObject,'TraceMask');
     if ~isempty(points)
      setappdata(hObject,'TracePoints',points(1:end-1,:));
      setappdata(hObject,'TracePath',path(1:end-1));
      setappdata(hObject,'TraceMask',mask(1:end-1));
      
      roi=getappdata(hObject,'uroi');
    %  roi=setv_points(roi,path{end},0);
      roi=setv_points(roi,mask{end},0);
      
      for i=1:length(path)-1
          roi=setv_points(roi,path{i},c.vpath);
      end
        
      
      setappdata(hObject,'uroi',roi);
      
      setappdata(handles.figure1,'isSaved',false);

      showImages(handles);
     end
     
 elseif eventdata.Character=='f' % finish the current path
     path=getappdata(hObject,'TracePath'); 
      roi=getappdata(hObject,'uroi');
      val=get(handles.ROIVal,'String');
      val=str2double(val);
    for i=1:length(path)
        roi=setv_points(roi,path{i},val); 
    end
     setappdata(hObject,'TracePoints',[]);
     setappdata(hObject,'TracePath',[]); 
     setappdata(hObject,'TraceMask',[]); 
    
    if length(path)>0
      setappdata(hObject,'uroi',roi);
      
      setappdata(handles.figure1,'isSaved',false);

      showImages(handles);
    end
 elseif  length(eventdata.Modifier )==1 && eventdata.Key=='s'   && strcmp(eventdata.Modifier{1},'control')
     %saveROI_Callback(hObject, eventdata, handles,1);
     
 elseif length(eventdata.Modifier )==0 && eventdata.Character=='s' 
         
     isshow=get(handles.showROI,'Value');
     set(handles.showROI,'Value',~isshow);
     showImages(handles);
 

 elseif eventdata.Character=='o'
    LoadTrace_Callback(hObject, eventdata, handles);
 elseif eventdata.Character=='e'
   removeROI_Callback(hObject, eventdata, handles);
    
 end
 
 
