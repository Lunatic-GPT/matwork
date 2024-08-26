 function wheelScroll_Callback(hObject, eventdata,h)
 

handles=getappdata(hObject,'handles');
data=getappdata(hObject,'udata');

odata_rs=getappdata(hObject,'odata_rs');
ind_u2o=getappdata(hObject,'ind_u2o');

if size(data,3)>1  || size(odata_rs,3)>1
           layout=str2num(get(handles.Layout,'String'));

           maxNumCols=layout(2);
           maxNumRows=layout(1);

        str=get(handles.first_slice,'String');
        first_slice=str2num(str);

        if eventdata.VerticalScrollCount>0
            first_slice=first_slice+1;
        else
            first_slice=first_slice-1;
        end

        if first_slice< min(1,1-ind_u2o)
            first_slice=min(1,1-ind_u2o);
        end
        
        if first_slice> max(size(data,3),size(odata_rs,3)-ind_u2o)
            first_slice=max(size(data,3),size(odata_rs,3)-ind_u2o);
        end
        
        set(handles.first_slice,'String',first_slice);

elseif size(data,4)>1
    
    str=get(handles.cur_brik,'String');
    cur_brik=str2num(str);
    
    if eventdata.VerticalScrollCount>0
         cur_brik=cur_brik+1;
    else
         cur_brik=cur_brik-1;
    end
   
    cur_brik=mod(cur_brik-1,size(data,4))+1;
    set(handles.cur_brik,'String',num2str(cur_brik));
    
else
    return;
end


 showImages(handles); 

   
        %     figure(gh);
             