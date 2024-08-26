function Draw1_Callback(hObject, eventdata, handles)
% hObject    handle to Draw1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.figure1,'WindowButtonDownFcn','');
shape=get(handles.roiShapes,'Val');
   set(gcf,'Pointer', 'Crosshair');
if shape==1
 
    h = imfreehand(gca);
   
elseif shape==2
    h = impoly(gca,[]);
elseif shape==3 || shape==4
    h = imellipse(gca);
    
    
    pos=getPosition(h);
    
    if pos(3)==0
        hroiSize=findobj(handles.figure1,'Tag','ROISize');
        roiSize=get(hroiSize,'String');
        
        pos(3:4)=str2num(roiSize);
        pos(1:2)=pos(1:2)-pos(3:4)/2-0.5;
    end
    if shape==3
        pos(3:4)=min(pos(3:4));
        setPosition(h,pos);
        
        setFixedAspectRatioMode(h,true);
    end
elseif shape==5
    h = imrect(gca);
elseif shape==6
    h=drawpolyline();
 
end
 set(handles.figure1,'Pointer', 'Arrow');

% set(handles.figure1,'WindowButtonDownFcn',@MouseButtonDown_Callback);

if isempty(h)
    return;
end

slice=getappdata(gca,'slice');

i=1;
while ~isempty(findobj(gcf,'Tag',roi_pos_name(slice,i,shape)))
    i=i+1;
end

set(h,'Tag',roi_pos_name(slice,i,shape));


%    roi=saveROI_Callback(hObject, eventdata, handles,false);
%  p=fileparts(mfilename('fullpath'));
%   save(fullfile(p,'roi_autosave'),'roi');