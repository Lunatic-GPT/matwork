function [roi,roi_draw]=saveROI_Callback(hObject, eventdata, handles,dosave)
% hObject    handle to saveROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist('dosave','var')
    dosave=true;
end



voxsize=str2num(get(handles.u_voxsize,'String'));
center=str2num(get(handles.u_center,'String'));

ROIMatchOverlay=get(handles.ROIMatchOverlay,'Value');

if ROIMatchOverlay
    d=getappdata(handles.figure1,'odata');
    roi=getappdata(handles.figure1,'oroi');
    
    roifile=getappdata(handles.figure1,'oroifile');
    ufilename=getappdata(handles.figure1,'ofile');
    
else
    roi=getappdata(handles.figure1,'uroi');
    d=getappdata(handles.figure1,'udata');
    
    roifile=getappdata(handles.figure1,'uroifile');
    ufilename=getappdata(handles.figure1,'ufile');
end


sz=size(d);
if length(sz)==2
    sz(3)=1;
end

if isempty(roi)
    roi=zeros(sz(1:3));
end

applyTo = get(handles.listbox7,'Value');

if applyTo==1
    roiv=get(handles.ROIVal,'String');
    roiv=str2num(roiv);
else
    roiv=1;
end

roi_draw=roi*0;

if length(roiv)==1  % apply new values to roi
    slice=getappdata(gca,'slice');
    tmp=draw2roi(handles,sz,slice);
    roitmp=roi(:,:,slice);
    roitmp(tmp>0)=roiv;
    roi_draw(:,:,slice)=tmp;
    roi(:,:,slice)=roitmp;
else
    roi(roi==roiv(1))=roiv(2);
end


setappdata(handles.figure1,'isSaved',false);
%%
if dosave
    
    roi=inverseImageTransform(roi,handles);
  
    
    save_as=false;
    if ~strcmp(get(hObject,'Type'),'figure') && strcmp(get(hObject,'String'),'Save as')
        save_as=true;
    end
    
    [~,usuffix]=strtok2(ufilename,'.');
        
        
    if isempty(roifile)
        [fname,dname]=uiputfile(fullfile(pwd,['*',usuffix]),'ROI name');
        roifile=fullfile(dname,fname);
    elseif save_as
        [dname,~,ext]=fileparts(roifile);
        if isempty(dname)
            dname=pwd;
        end
        [fname,dname]=uiputfile(fullfile(dname,['*',ext]),'ROI name');
        roifile=fullfile(dname,fname);
    end
    

    if ~isempty(roifile)
        if strcmp(roifile(end-3:end),'.mat')
            d=roi;
            save(roifile,'d','voxsize','center');
            
        elseif strcmp(roifile(end-3:end),'.nii') || strcmp(roifile(end-2:end),'.gz')
            nii=getappdata( handles.figure1,'nii');
            if isempty(nii)
                warning('No nifti header found; no file saved');
                return;
            end
            
            nii.img=roi;
            
            nii.hdr.dime.dim(5)=size(nii.img,4);
            nii.hdr.dime.dim(1)=ndims(nii.img);

            if strcmp(roifile(end-2:end),'.gz')
                save_untouch_niigz(nii,roifile(1:end-7));
            else
                save_untouch_nii(nii,roifile(1:end-4));
            end
        else 
            roifile='';  %not valid file; 
        end
    end
    
    if ~isempty(roifile)
     setappdata(handles.figure1,'isSaved',true);
    end
    
    if ROIMatchOverlay
        setappdata(handles.figure1,'oroifile',roifile);
        
        
    else
        setappdata(handles.figure1,'uroifile',roifile);
    end
    showImages(handles);
end

