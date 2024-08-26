function varargout = imageGallery(varargin)
% IMAGEGALLERY M-file for imageGallery.fig
%      IMAGEGALLERY, by itself, creates a new IMAGEGALLERY or raises the existing
%      singleton*.
%
%      H = IMAGEGALLERY returns the handle to a new IMAGEGALLERY or the handle to
%      the existing singleton*.
% 
% 
% This GUI allows you to open several images and batch process all of them
% 



gui_Singleton = 1;

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageGallery_OpeningFcn, ...
                   'gui_OutputFcn',  @imageGallery_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT

% --- Executes just before imageGallery is made visible.


 function pos = myupdatefcn(obj,event_obj)
    % Display 'Time' and 'Amplitude'
    pos = get(event_obj,'Position');
    plot_TimeCourse;
    
    
    handles=getappdata(gcf,'handles');
    ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');
    
    if ROIMatchOverlay
        roi=getappdata(gcf,'oroi');  
    else
        roi=getappdata(gcf,'uroi');
    end
        
    udata=getappdata(gcf,'udata');
    
    odata=getappdata(gcf,'odata');
    
    slice=getappdata(gca,'slice');
    brik=get(handles.cur_brik,'String');
    brik=str2num(brik);
    
    obrik=get(handles.obrik,'String');
    obrik=str2num(obrik);
    
    mbrik=get(handles.roi_brik,'String');
    mbrik=str2num(mbrik);
    
    crop_ud=get(handles.cropud,'String');
    crop_ud=str2num(crop_ud);
    
    crop_lr=get(handles.croplr,'String');
    crop_lr=str2num(crop_lr);
    
    if slice>=1 && slice <=size(udata,3)
        str=sprintf('%6.5f',udata(pos(1)+crop_lr(1),pos(2)+crop_ud(1),slice,brik));
    else
        str='';
    end
    
    
    
    if ROIMatchOverlay
        ind_u2o=getappdata(gcf,'ind_u2o');
    else
        ind_u2o=0;
    end
    
    
    if ~isempty(odata)
      str=sprintf('%s; o=%5.4f',str,odata(pos(1)+crop_lr(1),pos(2)+crop_ud(1),slice+ind_u2o,obrik));
    end
    
    if ~isempty(roi)
      str=sprintf('%s; roi=%d',str,roi(pos(1)+crop_lr(1),pos(2)+crop_ud(1),slice+ind_u2o,mbrik));  
    end
   
    pos=sprintf('%d %d (%s)',pos(1)+crop_lr(1),pos(2)+crop_ud(1),str);
     
  %  txt = {['Time: ',num2str(pos(1))],['Amplitude: ',num2str(pos(2))]};
 

function imageGallery_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for imageGallery

handles.output = hObject;
setappdata(handles.figure1,'actionItems',[handles.imageSelection]);
dcm=datacursormode(handles.figure1);
set(handles.figure1,'toolbar','figure');
set(handles.figure1,'WindowScrollWheelFcn',@wheelScroll_Callback);

%set(handles.figure1,'WindowButtonDownFcn',@MouseButtonDown_Callback);

setappdata(handles.figure1,'handles',handles);

set(dcm,'Enable','on');
set(dcm,'UpdateFcn',@myupdatefcn);

% Update handles structureuser
guidata(hObject, handles);
% UIWAIT makes imageGallery wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imageGallery_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in imageSelection.
function imageSelection_Callback(hObject, eventdata, handles)

    c_dir= pwd;

   % if exist('img_history.mat','file')
   
   mname=mfilename('fullpath');
  dname=fileparts(mname);
  
   if exist(fullfile(dname,'imageSelectionExts.mat'),'file')
       load(fullfile(dname,'imageSelectionExts.mat'),'exts');
   else
       exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
   end
   
  [files, dir_Name] = uigetfile(exts,'Load Underlay','MultiSelect','on');
   if ~dir_Name
        return;
    end
  [tmp,tmp2,ext]=fileparts(files);
  
  exts=unique(cat(1,{['*',ext]},exts),'stable');
  
  save(fullfile(dname,'imageSelectionExts.mat'),'exts');
  
  cd(c_dir);
  
  
  
    if ~iscell(files)
        temp = files;
        files = cell(1);
        files{1} = temp;
    end
    
    if length(files)>1
        error('select only one file');
    end
   
    
    setappdata(handles.figure1,'ufile',files{1});
    setappdata(handles.figure1,'dname',dir_Name);
    %setappdata(handles.figure1,
    
    if strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
        files_tmp=dir([dir_Name,'/*',ext]);
        if length(files_tmp)>1
            ns=inputdlg(sprintf('Input Number of slices (1-%d)',length(files_tmp)));
        else
            ns={'1'};
        end
        d=ri(fullfile(dir_Name,files{1}),str2num(ns{1}));
    elseif strcmp(ext,'.mat')
        try
            voxsize= ri(fullfile(dir_Name,files{1}),[],[],'voxsize');
            set(handles.u_voxsize,'String',num2str(voxsize));
            center= ri(fullfile(dir_Name,files{1}),[],[],'center');
            set(handles.u_center,'String',num2str(center));
            
            d= ri(fullfile(dir_Name,files{1}),[],[],'d');
            
        catch
            d=ri(fullfile(dir_Name,files{1}));
        end
        
    else
        
        d=ri(fullfile(dir_Name,files{1}));
    end
    
    d=double(d);
    setappdata(handles.figure1,'udata',d);
    setappdata(handles.figure1,'odata_rs',[]); % overlay data may need to be resampled since underlay is changed. 
    setappdata(handles.figure1,'odata',[]);
    
    if isappdata(handles.figure1,'udata2')
       rmappdata(handles.figure1,'udata2');
    end
    
    set(handles.showFull,'Value',true);
    
    dim=showImages(handles);
    
    set(handles.showFull,'Value',false);
    
    if length(dim)<4;
        dim(4)=1;
    end
  set(handles.ts_brik,'String',sprintf('1 %d',dim(4)));
  cbrik=get(handles.cur_brik,'String');
  if str2num(cbrik)>dim(4)
      set(handles.cur_brik,'String',num2str(dim(4)));
  end
  
  if  strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
     [vox,center]=dcmDimCenter(dir_Name);
    
     vox_str=sprintf('%5.3f,',vox);
     set(handles.u_voxsize,'String',vox_str(1:end-1));
     center_str=sprintf('%5.3f,',center);
    set(handles.u_center,'String',center_str(1:end-1));
      
  end

% --- Executes on button press in closeAll.
function closeAll_Callback(hObject, eventdata, handles)
%     imtool close all
    h = findobj(0,'type','figure');
    f = find (h ==handles.figure1);
    close(h(f));
% hObject    handle to closeAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
    showImages(handles);
%     hAx = getappdata(handles.figure1,'hAxes');
%     rows = getappdata(handles.figure1,'rows');
%     cols = getappdata(handles.figure1,'cols');
%     rPitch = getappdata(handles.figure1,'rPitch');
%     scale = 1-getappdata(handles.figure1,'panelMargin');
%     
%        [c r] = ind2sub([maxNumCols maxNumRows],ind);
%         x = (c-1)*cPitch;
%         y = 1-pMargin-r*rPitch;
%         
%     for i=1:length(hAx)
%         [c r] = ind2sub([cols rows],i);
%         
%         pos = get(hAx(i),'Position');
%         val = get(hObject,'value');
%        % val=round(val);
%        % set(hObject,'value',val);
%         
%         top = rPitch*rows-val * (rPitch*rows-scale);
%         pos(2) = top-r*rPitch;
%         set(hAx(i),'position',pos);
%     end
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_roi.
function loadroi_Callback(hObject, eventdata, handles)
% hObject    handle to load_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   mname=mfilename('fullpath');
  dname_img=fileparts(mname);
  
   if exist(fullfile(dname_img,'roiExts.mat'),'file')
       load(fullfile(dname_img,'roiExts.mat'),'exts');
   else
       exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
   end
   
   [fname,dname]=uigetfile(exts,'Select the roi file','ROIs.mat');

   if ~dname
        return;
   end
    
  [tmp,tmp2,ext]=fileparts(fname);
  
  exts=unique(cat(1,{['*',ext]},exts),'stable');
  
  save(fullfile(dname_img,'roiExts.mat'),'exts');
  
  roiMatchOverlay = get(handles.ROIMatchOverlay,'Value');
  
  u_voxsize=str2num(get(handles.u_voxsize,'String'));
  u_center=str2num(get(handles.u_center,'String'));
  
  ors_center=u_center;
  ors_center(3)=get_odata_rs_slice_center(handles);
  
  odata_rs=getappdata(gcf,'odata_rs');
  udata=getappdata(gcf,'udata');
  
  
  if strcmp(ext,'.mat')
      try
          voxsize= ri(fullfile(dname,fname),[],[],'voxsize');
          center= ri(fullfile(dname,fname),[],[],'center');
          roi = ri(fullfile(dname,fname),[],[],'d');
      catch
          roi =ri(fullfile(dname,fname));
          voxsize=u_voxsize;
          if ~roiMatchOverlay
              center = u_center;                       
          else
              center = ors_center;
          end        
      end
      
  else
      roi=ri(fullfile(dname,fname));

          voxsize=u_voxsize;
          if ~roiMatchOverlay
              center = u_center;                       
          else
              center = ors_center;
          end        

      
  end
  
  if roiMatchOverlay
      
      if size(roi,3)~=size(odata_rs,3)  || abs(ors_center(3)-center(3))>0.1
          error('This is not supported yet');
      end
      
      if (ndims(odata_rs(:,:,:,1))~=ndims(roi(:,:,:,1)) ...
              || any(size(odata_rs(:,:,:,1))~=size(roi(:,:,:,1)))) ...
              || any(abs(u_voxsize-voxsize)>0.01) || any(ors_center(1:2)~=center(1:2))
          [roi]=reslice_extend(roi,voxsize,u_voxsize,size(odata_rs(:,:,:,1)),ors_center-center);  %% odata_rs also centered on u_center
      end
      setappdata(handles.figure1,'oroi',roi);
      setappdata(handles.figure1,'oroifile',fname);
  else
      if size(roi,3)~=size(udata,3)  ||  abs(u_center(3)-center(3))>0.1
          error('This is not supported yet');
      end
      
      if (ndims(udata(:,:,:,1))~=ndims(roi(:,:,:,1)) ...
              || any(size(udata(:,:,:,1))~=size(roi(:,:,:,1)))) ...
              || any(abs(u_voxsize-voxsize)>0.01) || any(abs(u_center(1:2)-center(1:2))>0.01)
          [roi]=reslice_extend(roi,voxsize,u_voxsize,size(udata(:,:,:,1)),u_center-center);    
      end
      
      setappdata(handles.figure1,'uroifile',fname);
      setappdata(handles.figure1,'uroi',roi);
  end
  

showImages(handles);

set(handles.clearROI,'Enable','on');



function  add_contour(varargin)
        

hAx = varargin{1};
ch=findobj(hAx,'Type','image');


img = get(ch,'CData');
%img = floor(img*100/max(img(:)));
img(img<=2) = 2;
ca=caxis(hAx);
cm=gray(ca(2));
cm(1,:)=[0,1,0];
colormap(hAx,cm);
for i=2:nargin
     
     gap_BW=bwmorph(varargin{i},'remove');
     img(gap_BW>0) = 0;     
     
end

set(ch,'CData',img);



% --- Executes on button press in copy_zoom.
function copy_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to copy_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xlm=get(gca,'XLim');
ylm=get(gca,'YLim');


   hAxes = getappdata(handles.figure1,'hAxes');
    for i=1:length(hAxes)
    xlim(hAxes(i),xlm);
    ylim(hAxes(i),ylm);

    end


% --- Executes on button press in zoom_in.
function zoom_in_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=get(hObject,'Value');
h= zoom;
set(handles.pan,'Value',0);
set(handles.zoom_out,'Value',0);
if v==1
 set(h,'Direction','in','Enable','on');
else
 set(h,'Enable','off');
end

% --- Executes on button press in zoom_out.
function zoom_out_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=get(hObject,'Value');
h= zoom;
set(handles.pan,'Value',0);
set(handles.zoom_in,'Value',0);
if v==1
 set(h,'Direction','out','Enable','on');
else
 set(h,'Enable','off');
end


% --- Executes on button press in pan.
function pan_Callback(hObject, eventdata, handles)
% hObject    handle to pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v=get(hObject,'Value');
set(handles.zoom_out,'Value',0);
set(handles.zoom_in,'Value',0);
if v==1
    pan on;
else
    pan off;
end




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
showImages(handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
showImages(handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cur_brik_Callback(hObject, eventdata, handles)
% hObject    handle to cur_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cur_brik as text
%        str2double(get(hObject,'String')) returns contents of cur_brik as a double
showImages(handles);

% --- Executes during object creation, after setting all properties.
function cur_brik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cur_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ts_brik_Callback(hObject, eventdata, handles)
% hObject    handle to ts_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_brik as text
%        str2double(get(hObject,'String')) returns contents of ts_brik as a double

    showImages(handles);

% --- Executes during object creation, after setting all properties.
function ts_brik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

showImages(handles);


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



 c_dir= pwd;

   mname=mfilename('fullpath');
  dname_img=fileparts(mname);
  
   if exist(fullfile(dname_img,'overlayExts.mat'),'file')
       load(fullfile(dname_img,'overlayExts.mat'),'exts');
   else
       exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
   end
   
   [files, dir_Name]=uigetfile(exts,'Select the overlay file','ROIs.mat');

   if ~dir_Name
        return;
   end
    
  [tmp,tmp2,ext]=fileparts(files);
  
  exts=unique(cat(1,{['*',ext]},exts),'stable');
  
  save(fullfile(dname_img,'overlayExts.mat'),'exts');
  
    
  cd(c_dir);
  
 if iscell(files)
     error('select only one file');
 end
 
  if ~dir_Name
        return;
  end
    
  if  strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
      [vox,center]=dcmDimCenter(dir_Name);
      
      vox_str=sprintf('%5.3f,',vox);
      set(handles.o_voxsize,'String',vox_str(1:end-1));
      center_str=sprintf('%5.3f,',center);
      set(handles.o_center,'String',center_str(1:end-1));
      
  end
  
    setappdata(handles.figure1,'ofile',files);
    setappdata(handles.figure1,'odir_Name',dir_Name);
 
    if strcmp(ext,'.mat')
        try
            voxsize= ri(fullfile(dir_Name,files),[],[],'voxsize');
            set(handles.o_voxsize,'String',num2str(voxsize));
            center= ri(fullfile(dir_Name,files),[],[],'center');
            set(handles.o_center,'String',num2str(center));
            
            d= ri(fullfile(dir_Name,files),[],[],'d');
            
        catch
            d=ri(fullfile(dir_Name,files));
        end
        
    else
        d=read_afni_sdt_images(fullfile(dir_Name,files));
    end
    
    d=double(d);
    
    rng=min_max(d);
    set(handles.edit12,'String',num2str(rng));
    setappdata(handles.figure1,'odata',d);
    setappdata(handles.figure1,'odata_rs',[]);  
    showImages(handles);
%    set(handles.copy_zoom,'Enable','on');
  
    
  
  uf=getappdata(handles.figure1,'ufile');
  
  if isappdata(handles.figure1,'mfiles');
      mfile=getappdata(handles.figure1,'mfiles');
      set(handles.figure1,'Name',sprintf('%s/%s/%s',uf,files,mfile));
  else   
   set(handles.figure1,'Name',sprintf('%s/%s',uf,files));
  end







function obrik_Callback(hObject, eventdata, handles)
% hObject    handle to obrik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of obrik as text
%        str2double(get(hObject,'String')) returns contents of obrik as a double


    showImages(handles);
% --- Executes during object creation, after setting all properties.
function obrik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obrik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


 c_dir= pwd;

    
  [files, dir_Name] = uigetfile({'data-set_*.xml';'*.hdr';'*.mat';'*.HEAD';'*.sdt'},'MultiSelect','off');
  cd(c_dir);
  
    if ~dir_Name
        return;
    end
    
    setappdata(handles.figure1,'mfile',files);
    setappdata(handles.figure1,'mdir_Name',dir_Name);
   
    
    showImages(handles);
%    set(handles.copy_zoom,'Enable','on');
  
    
  
  uf=getappdata(handles.figure1,'ufile');
  
  of=getappdata(handles.figure1,'ofile');
  set(handles.figure1,'Name',sprintf('%s/%s/%s',uf,of,files));
  set(handles.checkbox2,'Value',0); 

function mrange_Callback(hObject, eventdata, handles)
% hObject    handle to mrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mrange as text
%        str2double(get(hObject,'String')) returns contents of mrange as a double

    showImages(handles);

% --- Executes during object creation, after setting all properties.
function mrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
 if get(handles.showOverlay,'Value')==1
    showImages(handles);
 end

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


    showImages(handles);
% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3

    showImages(handles);

% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

showImages(handles);
% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showOverlay.
function showOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to showOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showOverlay
   % if get(handles.showOverlay,'Value')==1
    showImages(handles);
    %end
    
    
    



function nvox_thr_Callback(hObject, eventdata, handles)
% hObject    handle to nvox_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nvox_thr as text
%        str2double(get(hObject,'String')) returns contents of nvox_thr as a double

    showImages(handles);
    
% --- Executes during object creation, after setting all properties.
function nvox_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nvox_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function croplr_Callback(hObject, eventdata, handles)
% hObject    handle to croplr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of croplr as text
%        str2double(get(hObject,'String')) returns contents of croplr as a double
    showImages(handles);

% --- Executes during object creation, after setting all properties.
function croplr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to croplr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cropud_Callback(hObject, eventdata, handles)
% hObject    handle to cropud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cropud as text
%        str2double(get(hObject,'String')) returns contents of cropud as a double
    showImages(handles);

% --- Executes during object creation, after setting all properties.
function cropud_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 hAxes = getappdata(handles.figure1,'hAxes');

   ssl=get(handles.saveslices,'String');

   ssl=str2num(ssl);
   
   saveds=[];
for ia=1:length(hAxes)
    c_image=getappdata(hAxes(ia),'image');

    nd=get(handles.magn,'String');
    nd=str2double(nd);
    cm=getappdata(hAxes(ia),'colormap');
    slice=getappdata(hAxes(ia),'slice');
    if ~isempty(slice)&&any(slice==ssl)
        saveds(end+1)=slice;
    else
        continue;
    end
    
    [c_image,cm]=add_roi_contour(c_image,cm,hAxes(ia));
    c_image(isnan(c_image))=0;    
    y=zeros(nd,size(c_image,1),nd,size(c_image,2),3);
    for i=1:size(c_image,1)
        for j=1:size(c_image,2)
            ind=round(c_image(i,j));
            if ind<1
                ind=1;
            elseif ind>size(cm,1)
                ind=size(cm,1);
            else
            end
            tmp=shiftdim(cm(ind,:),-3);
         y(:,i,:,j,:)=repmat(tmp,[nd,1,nd,1,1]);
        end
    end
    y=reshape(y,[size(c_image,1)*nd,size(c_image,2)*nd,3]);
  %  fname=sprintf('%s_%d_%d_%d_%dX%d.tif',strtok(get(handles.edit16,'String'),'.'),ud(1),ud(2),lr(1),lr(2),slice);
    fname=sprintf('%s_X%d.tif',strtok(get(handles.edit16,'String'),'.'),slice);
    imwrite(y,fname,'TIFF','Resolution',300);
    fprintf('%s saved\n',fname);
end
set(handles.saveslices,'String',num2str(saveds));
%set(hObject,'Enable','on');


function [img,cm]= add_roi_contour(img,cm,hAx)
        
    
    slice=getappdata(hAx,'slice');
    img(img>size(cm,1))=size(cm,1);
    
    roi=saveROI_Callback(hObject, eventdata, handles,false);
     
    gap_BW=bwmorph(roi(:,:,slice),'remove');
      
    img(gap_BW>0)=size(cm,1)+1;
    
    %cm=cat(1,cm,[1,1,1]);
    if sum(gap_BW(:))>0
        cm=cat(1,cm,[0,1,0]);
    end





function saveslices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in clearROI.
function clearROI_Callback(hObject, eventdata, handles)
% hObject    handle to clearROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


setappdata(handles.figure1,'uroi',[]);

setappdata(handles.figure1,'oroi',[]);

udata=getappdata(handles.figure1,'udata');
ind_u2o=getappdata(handles.figure1,'ind_u2o');
odata=getappdata(handles.figure1,'odata');

smin=1+min([0,ind_u2o]);
smax=max([size(odata,3)-ind_u2o,size(udata,3)]);

for i=smin:smax
    
    
    j=1;
    while isappdata(handles.figure1,roi_pos_name(i,j))
      rmappdata(handles.figure1,roi_pos_name(i,j));
      j=j+1;
    end
    
    for k=1:4
        j=1;
        while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(i,j,k)))
            delete(findobj(handles.figure1,'Tag',roi_pos_name(i,j,k)));
            j=j+1;
        end
    end
    
    
end

showImages(handles);

        

function Draw1_Callback(hObject, eventdata, handles)
    % hObject    handle to Draw1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
   % set(handles.figure1,'WindowButtonDownFcn','');
  shape=get(handles.roiShape,'Val');
  
  if shape==1
      h = impoly(gca,[]);
  elseif shape==2 || shape==3
      h = imellipse(gca);
      pos=getPosition(h);
      
      if shape==2
          pos(3:4)=min(pos(3:4));
          setPosition(h,pos);
          
          setFixedAspectRatioMode(h,true);
      end
  elseif shape==4
      h = imrect(gca);
      
  end
  
    
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
    roi=saveROI_Callback(hObject, eventdata, handles,false);
    p=fileparts(mfilename('fullpath'));
    save(fullfile(p,'roi_autosave'),'roi');
    
    
 % setappdata(gcf,roi_pos_name(slice,i),roi_pos);
% --- Executes on button press in saveROI.
function roi=saveROI_Callback(hObject, eventdata, handles,dosave)
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
    ind_u2o=getappdata(gcf,'ind_u2o');
    center(3)=get_odata_rs_slice_center(handles);
else
    roi=getappdata(handles.figure1,'uroi');
    d=getappdata(handles.figure1,'udata');
end

sz=size(d);
if length(sz)==2
    sz(3)=1;
end
              
if isempty(roi)
    roi=zeros(sz(1:3));
end

roiv=get(handles.ROIVal,'String');
roiv=str2num(roiv);
for i=1:size(roi,3)


    
    if ROIMatchOverlay
        slice=i-ind_u2o;
    else
        slice=i;
    end
  
    for k=1:4
        j=1;
        while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(slice,j,k)))
            h=findobj(handles.figure1,'Tag',roi_pos_name(slice,j,k));
            api = iptgetapi(h);
            roi_pos = api.getPosition();
            ud=get(handles.cropud,'String');
            ud=str2num(ud);
            lr=get(handles.croplr,'String');
            lr=str2num(lr);
            
            roi_pos(:,1)=roi_pos(:,1)+lr(1);
            roi_pos(:,2)=roi_pos(:,2)+ud(1);
            
            if k==1
                tmp=roipoly(zeros(sz(1:2)),roi_pos(:,2),roi_pos(:,1));  %x, y
            elseif k==2
                
                [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % x, y
                orig=roi_pos(1:2)+roi_pos(3:4)/2;
                tmp=sqrt((xx-orig(2)).^2+(yy-orig(1)).^2)<roi_pos(3)/2;
            elseif k==3
                
                [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % x, y
                orig=roi_pos(1:2)+roi_pos(3:4)/2;
                tmp=sqrt((xx-orig(2)).^2/roi_pos(4)^2+(yy-orig(1)).^2/roi_pos(3)^2)<1/4;
                
            elseif k==4
               
               posx=[roi_pos(1),roi_pos(1)+roi_pos(3),roi_pos(1)+roi_pos(3),roi_pos(1)];
               posy=[roi_pos(2),roi_pos(2),roi_pos(2)+roi_pos(4),roi_pos(2)+roi_pos(4)];  
               tmp=roipoly(zeros(sz(1:2)),posx,posy);
            end
            roitmp=roi(:,:,i);
            roitmp(tmp>0)=roiv;         
            roi(:,:,i)=roitmp;
            j=j+1;      
        end
        
    end
  %{  
    j=1;
    while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(slice,j,true)))
        h=findobj(handles.figure1,'Tag',roi_pos_name(slice,j,true));
        api = iptgetapi(h);
        roi_pos = api.getPosition();
        ud=get(handles.cropud,'String');
        ud=str2num(ud);
        lr=get(handles.croplr,'String');
        lr=str2num(lr);
                 
        roi_pos(1)=roi_pos(1)+lr(1);
                 
        roi_pos(2)=roi_pos(2)+ud(1);
      
       
        roitmp=roi(:,:,i);
        
        [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % x, y
        orig=roi_pos(1:2)+roi_pos(3:4)/2;
        roitmp(sqrt((xx-orig(2)).^2+(yy-orig(1)).^2)<roi_pos(3)/2)=roiv;
        roi(:,:,i)=roitmp;
      j=j+1;
    end
    %}
end



    %%
if dosave
[fname,dname]=uiputfile('*.mat','ROI name');
d=roi;
save(fullfile(dname,fname),'d','voxsize','center');
end

function res=get_odata_rs_slice_center(handles)
    
    o_center=str2num(get(handles.o_center,'String'));    
    u_voxsize=str2num(get(handles.u_voxsize,'String'));
    u_center=str2num(get(handles.u_center,'String'));
    
    center_u2o=u_center-o_center;
    
    nvox2shft = (center_u2o./u_voxsize);
    
    
    res=o_center(3)+(round(nvox2shft(3))-nvox2shft(3)).*u_voxsize(3);
    


            
% --- Executes on button press in removeROI.
function removeROI_Callback(hObject, eventdata, handles)
% hObject    handle to removeROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%set(handles.figure1,'WindowButtonDownFcn','');

h = impoly(gca,[]);

%set(handles.figure1,'WindowButtonDownFcn',@MouseButtonDown_Callback);

if isempty(h)
    return;
end
api = iptgetapi(h);
api.setColor('green');
roi_pos = api.getPosition();
ud=get(handles.cropud,'String');
ud=str2num(ud);
lr=get(handles.croplr,'String');
lr=str2num(lr);

roi_pos(:,1)=roi_pos(:,1)+lr(1);

roi_pos(:,2)=roi_pos(:,2)+ud(1);

d=getappdata(handles.figure1,'udata');
sz=size(d);
slice=getappdata(gca,'slice');

ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');
ind_u2o=getappdata(handles.figure1,'ind_u2o');
m =roipos2roi(roi_pos,sz(1:2));

if ROIMatchOverlay

    slice=slice+ind_u2o;
    
end

if ROIMatchOverlay
 roi2=getappdata(handles.figure1,'oroi');
else
 roi2=getappdata(handles.figure1,'uroi');
end

roi2(:,:,slice)=double(roi2(:,:,slice)).*double(m==0);

if ROIMatchOverlay
    setappdata(handles.figure1,'oroi',roi2);
else
    setappdata(handles.figure1,'uroi',roi2);
end

showImages(handles);

  

function first_slice_Callback(hObject, eventdata, handles)
% hObject    handle to first_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_slice as text
%        str2double(get(hObject,'String')) returns contents of first_slice as a double


   
   layout=str2num(get(handles.Layout,'String'));
   
   maxNumCols=layout(2);
   maxNumRows=layout(1);
   
str=get(handles.first_slice,'String');
first_slice=str2num(str);

   showImages(handles); 
    
% --- Executes during object creation, after setting all properties.
function first_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function magn_Callback(hObject, eventdata, handles)
% hObject    handle to magn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magn as text
%        str2double(get(hObject,'String')) returns contents of magn as a double


% --- Executes during object creation, after setting all properties.
function magn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in detectoROI.
function detectoROI_Callback(hObject, eventdata, handles)
% hObject    handle to detectoROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ROIMatchOverlay,'Value',true);
calcROI_Callback(hObject, eventdata, handles);


function calcROI_Callback(hObject, eventdata, handles)
        
sig=get(handles.sigma,'String');

sig=str2num(sig);

%%
matchOverlay=get(handles.ROIMatchOverlay,'Value');

if matchOverlay
    d=getappdata(handles.figure1,'odata_rs');
    ind_u2o=getappdata(handles.figure1,'ind_u2o');
    brik=get(handles.obrik,'String');
    brik=str2num(brik);
else
    d=getappdata(handles.figure1,'udata');
    ind_u2o=0;
    brik=get(handles.cur_brik,'String');
    brik=str2num(brik);
end

if matchOverlay
   roi=getappdata(handles.figure1,'oroi');
else
   roi=getappdata(handles.figure1,'uroi');
end

if isempty(roi) || size(roi,1)~=size(d,1) || size(roi,2)~=size(d,2) || size(roi,3)~=size(d,3) 
    roi=0*d(:,:,:,1);
end

%%

roival=get(handles.ROIVal,'String');
roival=str2num(roival);


%roival=ceil(roival/2)*2;
%roic=(roi==roival-1);
slice=getappdata(gca,'slice');

h=findobj(handles.figure1,'Tag',sprintf('roi_pos%d_1',slice));
api = iptgetapi(h);
roi_pos = api.getPosition();
ud=get(handles.cropud,'String');
ud=str2num(ud);
lr=get(handles.croplr,'String');
lr=str2num(lr);

roi_pos(:,1)=roi_pos(:,1)+lr(1);

roi_pos(:,2)=roi_pos(:,2)+ud(1);

roi_tmp=roipoly(0*d(:,:,1,1),roi_pos(:,2),roi_pos(:,1));

roi(:,:,slice+ind_u2o)=roi(:,:,slice+ind_u2o).*double(roi_tmp==0);

tmp=d(:,:,slice+ind_u2o,brik);
    
y=sort(tmp(roi_tmp>0));

 if sig>0
    mn=mean(y(round(end/8):end));
    sd=std(y(round(end/8):end));
   mask=(tmp>mn+sig*sd)&roi_tmp>0; 
 else   
     mn=mean(y(1:round(end/8*7)));
    sd=std(y(1:round(end/8*7)));
   
   mask=tmp<mn+sig*sd&roi_tmp>0; 
 end
%%
mc=clusterize2(mask,2);
roi(:,:,slice+ind_u2o)=roi(:,:,slice+ind_u2o)+mc.*mask;

if matchOverlay
    setappdata(handles.figure1,'oroi',roi);
else
    setappdata(handles.figure1,'uroi',roi);
end

showImages(handles);





function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in roiright.
function roiShift_Callback(hObject, eventdata, handles)
% hObject    handle to roiright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  ROIMatchOverlay=get(handles.ROIMatchOverlay,'Value');

  tag=get(hObject,'Tag');
  
if ROIMatchOverlay
    m=getappdata(handles.figure1,'oroi');
else 
    m=getappdata(handles.figure1,'uroi');
end
    
step=get(handles.roistep,'String');
step=str2num(step);

  slice=getappdata(gca,'slice');
  val=get(handles.ROIVal,'String');
  val=str2num(val);
  
  if ROIMatchOverlay
      ind_u2o=getappdata(gcf,'ind_u2o');
  else
      ind_u2o=0;
  end
  
  roiShiftAll=get(handles.roiShiftAll,'Value');
  roiShiftAllSlices=get(handles.roiShiftAllSlices,'Value');
  
  if roiShiftAllSlices
      slice=1:size(m,3);
      ind_u2o=0;
  end
  
  roiShiftVisible=get(handles.shiftVisibleROIs,'Value');
  if roiShiftVisible
    xl=xlim;
    yl=ylim;
    
    xl=round(xl(1)):round(xl(2));
  
    yl=round(yl(1)):round(yl(2));
  
    set(handles.roiShiftAll,'Value',false);
    roiShiftAll=false;
    val=unique(vec(m(xl,yl,slice+ind_u2o)));
    val(val==0)=[];
  end
  
  if roiShiftAll
      
   m2=m(:,:,slice+ind_u2o);
  else
   m2=0*m(:,:,slice+ind_u2o);
   
   for i=1:length(val) 
       m2=m2+val(i)*(m(:,:,slice+ind_u2o)==val(i));
   end
   
  end
  
  m3=m(:,:,slice+ind_u2o);
  
  m3(m2>0)=0;
  
  if strcmp(tag,'roiright')
    m(:,:,slice+ind_u2o)=m3+circshift(m2,[step,0,0]);
  elseif strcmp(tag,'roileft')
      m(:,:,slice+ind_u2o)=m3+circshift(m2,[-step,0,0]);
  elseif strcmp(tag,'roiup')
      m(:,:,slice+ind_u2o)=m3+circshift(m2,[0,-step,0]);
  else
      m(:,:,slice+ind_u2o)=m3+circshift(m2,[0,step,0]);
  end
      
  
  
  ax=gca;
  if ROIMatchOverlay
      setappdata(handles.figure1,'oroi',m);
  else
      setappdata(handles.figure1,'uroi',m);
  end
      
  showImages(handles);
  axes(ax);
  
% --- Executes on button press in roi_fliplr.
function roi_fliplr_Callback(hObject, eventdata, handles)
% hObject    handle to roi_fliplr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ax=gca;
m=getappdata(handles.figure1,'roi');

slice=getappdata(gca,'slice');


roival=get(handles.FlipROIVal,'String');
roival=str2num(roival);

tmp=flipdim(m(:,:,slice,:),1);

tmp(tmp>0)=roival;

m=mask_on_briks(m,tmp,slice,roival);

setappdata(handles.figure1,'roi',m);

showImages(handles);
axes(ax);


function roistep_Callback(hObject, eventdata, handles)
% hObject    handle to roistep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roistep as text
%        str2double(get(hObject,'String')) returns contents of roistep as a double


% --- Executes during object creation, after setting all properties.
function roistep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roistep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cross.
function cross_Callback(hObject, eventdata, handles)
% hObject    handle to cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m=getappdata(handles.figure1,'roi');
step=get(handles.roistep,'String');
step=str2num(step);

  %slice=getappdata(gca,'slice');
  
  m=circshift(m,[0,0,step,0]);
  ax=gca;
  setappdata(handles.figure1,'roi',m);
  
  showImages(handles);
  axes(ax);
  
function nCol_Callback(hObject, eventdata, handles)
% hObject    handle to nCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nCol as text
%        str2double(get(hObject,'String')) returns contents of nCol as a double
showImages(handles);

% --- Executes during object creation, after setting all properties.
function nCol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Layout_Callback(hObject, eventdata, handles)
% hObject    handle to Layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Layout as text
%        str2double(get(hObject,'String')) returns contents of Layout as a double
showImages(handles);


% --- Executes during object creation, after setting all properties.
function Layout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roistat.
function roistat_Callback(hObject, eventdata, handles)
% hObject    handle to roistat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    roi=saveROI_Callback(hObject, eventdata, handles,false);

    udata=getappdata(handles.figure1,'udata');

    ind=get(handles.cur_brik,'String');
    ind=str2num(ind);
    tmp=udata(:,:,:,ind);
    
    roiVal= get(handles.ROIVal,'String');
    roiVal=str2double(roiVal);
    ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');
    
    if ROIMatchOverlay
        ind_u2o=getappdata(gcf,'ind_u2o');    
        roiu=0*tmp;
        
        for i=1:size(roiu,3)
            
            if i+ind_u2o>1 && i+ind_u2o<=size(roi,3)
            roiu(:,:,i)=roi(:,:,i+ind_u2o);
            end
            
        end
    else
       roiu=roi;
    end
    
    if roiVal ==0
        mn=mean(tmp(roiu>0));
        sd=std(tmp(roiu>0));
        fprintf('Sum = %f; Mean = %f; Std = %f;\n', sum(tmp(roiu>0)),mn,sd);
        fprintf('Max/Min = %f/%f; %d voxels\n',max(tmp(roiu>0)),min(tmp(roiu>0)),sum(roiu(:)>0));
        
    else
        mn=mean(tmp(roiu==roiVal));
        sd=std(tmp(roiu==roiVal));
        % fprintf('Mean = %f; Std = %f; Max/Min = %f/%f; %d voxels\n',mn,sd,max(tmp(roi==roiVal)),min(tmp(roi==roiVal)),sum(roi(:)==roiVal));
        disp('----------------------------------------------------------------------------------------------');
        fprintf('%d voxels, Sum = %e; Mean = %e; Std = %e; Max/Min = %e/%e\n',sum(roiu(:)==roiVal), sum(tmp(roiu==roiVal)),mn,sd,...
            max(tmp(roiu==roiVal)),min(tmp(roiu==roiVal)));
        
        
    end
    
    
%msgbox(sprintf('Mean = %f\n std = %f',mn,sd));


%fprintf('Mean = \n');
%disp(mean_roi(udata,roi>0));


% --- Executes on button press in swapxy.
function swap_Callback(hObject, eventdata, handles)
% hObject    handle to swapxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag=get(hObject,'Tag');


if strcmp(tag,'swapxy')
    new_order=[2,1,3,4];
elseif strcmp(tag,'swapyz')
    new_order=[1,3,2,4];
elseif strcmp(tag,'swapxz')
    new_order=[3,2,1,4];
end
    
    
    
d=getappdata(handles.figure1,'udata');
d=permute(d,new_order);
setappdata(handles.figure1,'udata',d);

d=getappdata(handles.figure1,'odata');
d=permute(d,new_order);
setappdata(handles.figure1,'odata',d);

d=getappdata(handles.figure1,'umdata');
d=permute(d,new_order);
setappdata(handles.figure1,'umdata',d);


d=getappdata(handles.figure1,'odata_rs');
d=permute(d,new_order);
setappdata(handles.figure1,'odata_rs',d);

uroi=getappdata(handles.figure1,'uroi');
uroi=permute(uroi,new_order);
setappdata(handles.figure1,'uroi',uroi);

oroi=getappdata(handles.figure1,'oroi');
oroi=permute(oroi,new_order);
setappdata(handles.figure1,'oroi',oroi);


o_voxsize=str2num(get(handles.o_voxsize,'String'));
 set(handles.o_voxsize,'String',num2str(o_voxsize(new_order(1:3))));

 o_center=str2num(get(handles.o_center,'String'));
 set(handles.o_center,'String',num2str(o_center(new_order(1:3))));
            
 u_voxsize=str2num(get(handles.u_voxsize,'String'));
 set(handles.u_voxsize,'String',num2str(u_voxsize(new_order(1:3))));
 
 u_center=str2num(get(handles.u_center,'String'));
 set(handles.u_center,'String',num2str(u_center(new_order(1:3))));
 

  showImages(handles);




% --- Executes on button press in fliplr_image.
function flip_image_Callback(hObject, eventdata, handles)
% hObject    handle to fliplr_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tag=get(hObject,'Tag');
if strcmp(tag,'fliplr_image')
    dim=1;
else
    dim=2;
end

d=getappdata(handles.figure1,'udata');
d=flip(d,dim);
setappdata(handles.figure1,'udata',d);

d=getappdata(handles.figure1,'odata');
d=flip(d,dim);
setappdata(handles.figure1,'odata',d);

d=getappdata(handles.figure1,'umdata');
d=flip(d,dim);
setappdata(handles.figure1,'umdata',d);


d=getappdata(handles.figure1,'odata_rs');
d=flip(d,dim);
setappdata(handles.figure1,'odata_rs',d);

uroi=getappdata(handles.figure1,'uroi');
uroi=flip(uroi,dim);
setappdata(handles.figure1,'uroi',uroi);

oroi=getappdata(handles.figure1,'oroi');
oroi=flip(oroi,dim);
setappdata(handles.figure1,'oroi',oroi);

%need to recalculate voxsize and center for odata, but ok if not
%recalculate odata_rs.

 
showImages(handles);


function saveslices_Callback(hObject, eventdata, handles)
% hObject    handle to saveslices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveslices as text
%        str2double(get(hObject,'String')) returns contents of saveslices as a double


% --- Executes on button press in Go2Vox.
function Go2Vox_Callback(hObject, eventdata, handles)
% hObject    handle to Go2Vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    matchOverlay=get(handles.ROIMatchOverlay,'Value');

    if matchOverlay
    
    d=getappdata(handles.figure1,'oroi');

    else

    d=getappdata(handles.figure1,'uroi');
        
    end
val=get(handles.ROIVal,'String');
val=str2double(val);
os=get(handles.VoxOffset,'String');
os=str2double(os);
stp=get(handles.VoxStep,'String');
stp=str2double(stp);


ws=get(handles.WinSz,'String');
ws=str2double(ws);

ind=find(d==val);

sz=size(d);
if length(sz)==2
    sz(3)=1;
end
os=mod(os-1,length(ind))+1;

[i1,i2,i3]=ind2sub(sz,ind(os));

os=os+stp;
set(handles.VoxOffset,'String',num2str(os));

nlayout=get(handles.Layout,'String');
nlayout=str2num(nlayout);
nrow=nlayout(1);
ncol=nlayout(2);


set(handles.first_slice,'String',num2str(i3-floor((nrow*ncol)/2)));
first_slice_Callback(hObject, eventdata, handles);

showImages(handles);
xl=i1+[-1,1]/2*ws;
yl=i2+[-1,1]/2*ws;
if xl(1)<1
    xl(1)=1;
end
if yl(1)<1
    yl(1)=1;
end


xlim(xl);
ylim(yl);





function VoxOffset_Callback(hObject, eventdata, handles)
% hObject    handle to VoxOffsetText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VoxOffsetText as text
%        str2double(get(hObject,'String')) returns contents of VoxOffsetText as a double


% --- Executes during object creation, after setting all properties.
function VoxOffsetText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoxOffsetText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in showFull.
function showFull_Callback(hObject, eventdata, handles)
% hObject    handle to showFull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.showFull,'Value',true);
showImages(handles);
set(handles.showFull,'Value',false);
% Hint: get(hObject,'Value') returns toggle state of showFull



function WinSz_Callback(hObject, eventdata, handles)
% hObject    handle to WinSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WinSz as text
%        str2double(get(hObject,'String')) returns contents of WinSz as a double


% --- Executes during object creation, after setting all properties.
function WinSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WinSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VoxStep_Callback(hObject, eventdata, handles)
% hObject    handle to VoxStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VoxStep as text
%        str2double(get(hObject,'String')) returns contents of VoxStep as a double


% --- Executes during object creation, after setting all properties.
function VoxStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoxStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function VoxOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoxOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
if length(sz)==2;
    sz(3)=1;
end
if isempty(ind_u2o)
    ind_u2o=0;
end
if isempty(roi)
    
    roi=zeros(sz(1:3));
    
end


save ApplyROI_temp_before roi


roiv=get(handles.ROIVal,'String');
roiv=str2num(roiv);
                
for i=1:size(roi,3)
    j=1;
    while isappdata(handles.figure1,roi_pos_name(i-ind_u2o,j))
     roi_pos=getappdata(handles.figure1,roi_pos_name(i-ind_u2o,j));
      if ~isempty(roi_pos)
          
      roi(:,:,i)=roipoly(zeros(sz(1:2)),roi_pos(:,2),roi_pos(:,1)); 
      end
      j=j+1;
    end
    
    j=1;
    while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(i-ind_u2o,j)))
        h=findobj(handles.figure1,'Tag',roi_pos_name(i-ind_u2o,j));
        api = iptgetapi(h);
        roi_pos = api.getPosition();
        ud=get(handles.cropud,'String');
        ud=str2num(ud);
        lr=get(handles.croplr,'String');
        lr=str2num(lr);
                 
        roi_pos(:,1)=roi_pos(:,1)+lr(1);
                 
        roi_pos(:,2)=roi_pos(:,2)+ud(1);
                      
       
       tmp =roipos2roi(roi_pos,sz);
        roitmp=roi(:,:,i);
        roitmp(tmp>0)=roiv;
        
        roi(:,:,i)=roitmp;
        
        
      j=j+1;
      
      delete(h);
    end
    
    
    j=1;
    while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(i-ind_u2o,j,true)))
        h=findobj(handles.figure1,'Tag',roi_pos_name(i-ind_u2o,j,true));
        api = iptgetapi(h);
        roi_pos = api.getPosition();
        ud=get(handles.cropud,'String');
        ud=str2num(ud);
        lr=get(handles.croplr,'String');
        lr=str2num(lr);
                 
        roi_pos(1)=roi_pos(1)+lr(1);  %x
                 
        roi_pos(2)=roi_pos(2)+ud(1); %y
      
       
        roitmp=roi(:,:,i);
        
        [xx,yy]=meshgrid(1:size(roi,2),1:size(roi,1));  % y, x
        orig=roi_pos(1:2)+roi_pos(3:4)/2;
        
        roitmp(sqrt((xx-orig(2)).^2+(yy-orig(1)).^2)<roi_pos(3)/2)=roiv;
        
        roi(:,:,i)=roitmp;
      j=j+1;
      delete(h);
    end
    
    
end


if matchOverlay
 setappdata(handles.figure1,'oroi',roi);
else
 setappdata(handles.figure1,'uroi',roi);
end    
    
save ApplyROI_temp_after roi
showImages(handles);

function roi=roipos2roi(roi_pos,sz)
        
 tmp=roipoly(zeros(sz(1:2)),roi_pos(:,2),roi_pos(:,1));
        
      %  if sum(tmp(:))==0
            
            if size(roi_pos,1)==1
              tmp (round(roi_pos(1)),round(roi_pos(2)))=1;
            else
                for kk=1:size(roi_pos,1)-1
                    
                   n=2*sqrt((roi_pos(kk,1)-roi_pos(kk+1,1))^2+(roi_pos(kk,2)-roi_pos(kk+1,2))^2);
                   n=round(n);
                       ix=linspace(round(roi_pos(kk,1)),round(roi_pos(kk+1,1)),n);
                   
                       iy=linspace(round(roi_pos(kk,2)),round(roi_pos(kk+1,2)),n);
                   
                   for kkk=1:n
                     tmp (round(ix(kkk)),round(iy(kkk)))=1;
                   end
                end
            
            end
            
            
       % end
        
  roi=tmp;
        
function ROIVal_Callback(hObject, eventdata, handles)
% hObject    handle to ROIVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROIVal as text
%        str2double(get(hObject,'String')) returns contents of ROIVal as a double


% --- Executes during object creation, after setting all properties.
function ROIVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showROI.
function showROI_Callback(hObject, eventdata, handles)
% hObject    handle to showROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showROI
showImages(handles);



% --- Executes on button press in LineStat.
function LineStat_Callback(hObject, eventdata, handles)
% hObject    handle to LineStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gap=get(handles.roigap,'String');
gap=str2num(gap);

 a=getappdata(handles.figure1,'udata');
 m=getappdata(handles.figure1,'roi');
        
clr='rbg';
cnr=[];
figure;
n=max(m(:));

 for val=1:ceil(n/2)
        
      in=find(m==val*2-1 | m==val*2);
      if isempty(in)
          continue;
      end
    in2=find(m==val*2);
    
    in=sort(in);
    n0=find(in==in2);
    
    bl=[1:n0-gap,n0+gap:length(in)];
    
    inbl=in(bl);
    
       
    cnr(val)=(a(in2)-mean(a(inbl)))/std(a(inbl));
    
    
    plot(1:length(in),a(in),clr(val));
    hold on;
    plot(n0,a(in2),['o',clr(val)]);
    
 end
 
 cnr(cnr==0)=[];
 disp(mean(cnr));


function roigap_Callback(hObject, eventdata, handles)
% hObject    handle to roigap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roigap as text
%        str2double(get(hObject,'String')) returns contents of roigap as a double


% --- Executes during object creation, after setting all properties.
function roigap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roigap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roiHist.
function roiHist_Callback(hObject, eventdata, handles)
% hObject    handle to roiHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


roi=saveROI_Callback(hObject, eventdata, handles,false);

    udata=getappdata(handles.figure1,'udata');

    ind=get(handles.cur_brik,'String');
    ind=str2num(ind);
    tmp=udata(:,:,:,ind);
    
mn=mean(tmp(roi>0));
sd=std(tmp(roi>0));

%msgbox(sprintf('Mean = %f\n std = %f',mn,sd));

fprintf('Mean = %f; Std = %f; Max/Min = %f/%f; %d voxels\n',mn,sd,max(tmp(roi>0)),min(tmp(roi>0)),sum(roi(:)>0));
x=tmp(roi>0);

xl=min(x);
xu=max(x);

figure;hist(x,linspace(xl,xu,20));



% --- Executes on button press in ROIFlash.
function ROIFlash_Callback(hObject, eventdata, handles)
% hObject    handle to ROIFlash (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cValue=get(handles.showROI,'Value');
set(handles.showROI,'Value',~cValue);

showImages(handles);
pause(0.3);
set(handles.showROI,'Value',cValue);

showImages(handles);


% --- Executes on button press in rmCon.
function rmCon_Callback(hObject, eventdata, handles)
% hObject    handle to rmCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



roi=getappdata(handles.figure1,'roi');
 sz=size(roi);


h=impoly(gca);
seed=zeros(size(roi));

        api = iptgetapi(h);
        roi_pos = api.getPosition();
        ud=get(handles.cropud,'String');
        ud=str2num(ud);
        lr=get(handles.croplr,'String');
        lr=str2num(lr);
                 
        roi_pos(:,1)=roi_pos(:,1)+lr(1);
                 
        roi_pos(:,2)=roi_pos(:,2)+ud(1);
                  
        slice=getappdata(gca,'slice');
       seed(:,:,slice)=roipos2roi(roi_pos,sz(1:2));
       
        
        

      delete(h);
    

roi=remove_con(roi,seed);
setappdata(handles.figure1,'roi',roi);
showImages(handles);

function mask=remove_con(mask,seed)
         
    
    a=find(seed>0);
    
    for i=1:length(a)
        cluster_temp(i,:)=ind2subb(size(mask),a(i));
    end
    
    ivc=1;
    
    nvc=length(a);
    
         while ivc<=nvc 
        
          nb = find_neighbors(cluster_temp(ivc,:),mask);
        
          nnb = size(nb,1);   
         if nnb > 0 
          cluster_temp(nvc+1:nvc+nnb,:) = nb;
          nvc=nvc+nnb;
          for a=1:nnb
             mask(nb(a,1),nb(a,2),nb(a,3))=0;
          end
         end
        ivc = ivc+1;
         end
         
         
         
       
        
function nb = find_neighbors(pos,mask)

   sz = size(mask);
   if numel(sz) <3
       sz(3)=1;
   end
  nb=[];
  for i=-1:1
      for j=-1:1
          for k=-1:1
            
              pos2=pos+[i,j,k];
              if pos2(1)<1 || pos2(1)>sz(1) ||pos2(2)<1 || pos2(2)>sz(2)||pos2(3)<1 || pos2(3)>sz(3)
                  continue;
              end
              
              if mask(pos2(1),pos2(2),pos2(3))>0
                  nb(end+1,:)=pos2;
              end
          end
      end
  end
      


% --- Executes on button press in tFT.
function tFT_Callback(hObject, eventdata, handles)
% hObject    handle to tFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tFT


% --- Executes on button press in ROIFlash1.
function ROIFlash1_Callback(hObject, eventdata, handles)
% hObject    handle to ROIFlash1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

matchOverlay = get(handles.ROIMatchOverlay,'Value');
if matchOverlay
   roi=getappdata(handles.figure1,'oroi');
else
     roi=getappdata(handles.figure1,'uroi');
end

ROIVal=get(handles.ROIVal,'String');
ROIVal=str2num(ROIVal);

roi_orig=roi;

roi(roi==ROIVal)=0;

if matchOverlay
    setappdata(handles.figure1,'oroi',roi);
else
    setappdata(handles.figure1,'uroi',roi);
end

set(handles.showROI,'Value',true);

showImages(handles);

pause(0.3);

if matchOverlay
    setappdata(handles.figure1,'oroi',roi_orig);
else
    setappdata(handles.figure1,'uroi',roi_orig);
end
showImages(handles);







% --- Executes on button press in ROIClust.
function ROIClust_Callback(hObject, eventdata, handles)
% hObject    handle to ROIClust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');


if ROIMatchOverlay    
 roi=getappdata(handles.figure1,'oroi');
else
 roi=getappdata(handles.figure1,'uroi');
end
    
roi=clusterize2(roi>0,1,0);

if ROIMatchOverlay    
 setappdata(handles.figure1,'oroi',roi);
else
 setappdata(handles.figure1,'uroi',roi);
end


fprintf('Number of clusters = %d\n',max(roi(:)));
showImages(handles);

% --- Executes on button press in ROIMIP.
function ROIMIP_Callback(hObject, eventdata, handles)
% hObject    handle to ROIMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi=getappdata(handles.figure1,'roi');


val=get(handles.mipslice,'String');
val=str2num(val);

slice=getappdata(gca,'slice');
ind_u2o=getappdata(gcf,'ind_u2o');

roi2=roi;
roi2(:,:,slice+ind_u2o)=sum(roi(:,:,slice+val+ind_u2o),3);
setappdata(handles.figure1,'roi',roi2);

showImages(handles);



pause(0.3);

setappdata(handles.figure1,'roi',roi);

showImages(handles);





function mipslice_Callback(hObject, eventdata, handles)
% hObject    handle to mipslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mipslice as text
%        str2double(get(hObject,'String')) returns contents of mipslice as a double


% --- Executes during object creation, after setting all properties.
function mipslice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mipslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roiname_Callback(hObject, eventdata, handles)
% hObject    handle to roiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roiname as text
%        str2double(get(hObject,'String')) returns contents of roiname as a double


% --- Executes during object creation, after setting all properties.
function roiname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nroicolor_Callback(hObject, eventdata, handles)
% hObject    handle to nroicolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nroicolor as text
%        str2double(get(hObject,'String')) returns contents of nroicolor as a double

showImages(handles);
% --- Executes during object creation, after setting all properties.
function nroicolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nroicolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotts.
function plotts_Callback(hObject, eventdata, handles)
% hObject    handle to plotts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotts


% --- Executes on button press in SaveMosaic.
function SaveMosaic_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


 hAxes = getappdata(handles.figure1,'hAxes');

 Layout=get(handles.Layout,'String');
   layout=str2num(Layout);
   c_image=getappdata(hAxes(1),'image');
sz=size(c_image);
       nd=get(handles.magn,'String');
    nd=str2double(nd);
   yall=zeros(sz(1)*nd*layout(1),sz(2)*nd*layout(2),3);
   slice=[];
for ia=1:length(hAxes)
    c_image=getappdata(hAxes(ia),'image');


    cm=getappdata(hAxes(ia),'colormap');
    slice(end+1)=getappdata(hAxes(ia),'slice');
   
    
    [c_image,cm]=add_roi_contour(c_image,cm,hAxes(ia));
        
    y=zeros(nd,size(c_image,1),nd,size(c_image,2),3);
    for i=1:size(c_image,1)
        for j=1:size(c_image,2)
            ind=round(c_image(i,j));
            if isnan(ind)
                ind=1;
            end
            if ind<1
                ind=1;
            elseif ind>size(cm,1)
                ind=size(cm,1);
            else
            end
            tmp=shiftdim(cm(ind,:),-3);
         y(:,i,:,j,:)=repmat(tmp,[nd,1,nd,1,1]);
        end
    end
    y=reshape(y,[size(c_image,1)*nd,size(c_image,2)*nd,3]);
  %  fname=sprintf('%s_%d_%d_%d_%dX%d.tif',strtok(get(handles.edit16,'String'),'.'),ud(1),ud(2),lr(1),lr(2),slice);
 
  icol=mod(ia-1,layout(2))+1;
  irow=ceil(ia/layout(2));
  
  yall((irow-1)*sz(1)*nd+1:irow*nd*sz(1),(icol-1)*sz(2)*nd+1:icol*nd*sz(2),:)=y;
   
   
end

 fname=sprintf('%s_X%d_X%d.tif',strtok(get(handles.edit16,'String'),'.'),min(slice),max(slice));
    imwrite(yall,fname,'TIFF','Resolution',300);
    fprintf('%s saved\n',fname);


function odim_Callback(hObject, eventdata, handles)
% hObject    handle to odim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of odim as text
%        str2double(get(hObject,'String')) returns contents of odim as a double

setappdata(handles.figure1,'doresample',true);
showImages(handles);

% --- Executes during object creation, after setting all properties.
function odim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to odim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function centerShift_Callback(hObject, eventdata, handles)
% hObject    handle to centerShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of centerShift as text
%        str2double(get(hObject,'String')) returns contents of centerShift as a double

setappdata(handles.figure1,'doresample',true);
showImages(handles);
             
% --- Executes during object creation, after setting all properties.
function centerShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centerShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roiShiftAll.
function roiShiftAll_Callback(hObject, eventdata, handles)
% hObject    handle to roiShiftAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roiShiftAll


% --- Executes on button press in olup.
function olup_Callback(hObject, eventdata, handles)
% hObject    handle to olup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%{
step=get(handles.odataShiftStep,'String');
step=str2num(step);

%odata=getappdata(handles.figure1,'odata');
% odata=circshift(odata,-step,2);
%setappdata(handles.figure1,'odata',odata);

odata_rs=getappdata(handles.figure1,'odata_rs');
odata_rs=circshift(odata_rs,-step,2);
setappdata(handles.figure1,'odata_rs',odata_rs);

showImages(handles);
%}

% --- Executes on button press in oldown.
function oldown_Callback(hObject, eventdata, handles)
% hObject    handle to oldown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
step=get(handles.odataShiftStep,'String');
step=str2num(step);

%odata=getappdata(handles.figure1,'odata');
%odata=circshift(odata,step,2);
%setappdata(handles.figure1,'odata',odata);

odata_rs=getappdata(handles.figure1,'odata_rs');
odata_rs=circshift(odata_rs,step,2);
setappdata(handles.figure1,'odata_rs',odata_rs);

showImages(handles);


% --- Executes on button press in olleft.
function olleft_Callback(hObject, eventdata, handles)
% hObject    handle to olleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

step=get(handles.odataShiftStep,'String');
step=str2num(step);

%odata=getappdata(handles.figure1,'odata');
%odata=circshift(odata,-step,1);
%setappdata(handles.figure1,'odata',odata);

odata_rs=getappdata(handles.figure1,'odata_rs');
odata_rs=circshift(odata_rs,-step,1);
setappdata(handles.figure1,'odata_rs',odata_rs);

showImages(handles);

% --- Executes on button press in olright.
function olright_Callback(hObject, eventdata, handles)
% hObject    handle to olright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

step=get(handles.odataShiftStep,'String');
step=str2num(step);

%odata=getappdata(handles.figure1,'odata');
%odata=circshift(odata,step,1);
%setappdata(handles.figure1,'odata',odata);

odata_rs=getappdata(handles.figure1,'odata_rs');
odata_rs=circshift(odata_rs,step,1);
setappdata(handles.figure1,'odata_rs',odata_rs);

showImages(handles);


function roi_brik_Callback(hObject, eventdata, handles)
% hObject    handle to roi_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_brik as text
%        str2double(get(hObject,'String')) returns contents of roi_brik as a double

showImages(handles);

% --- Executes during object creation, after setting all properties.
function roi_brik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4
  showImages(handles);

% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5
showImages(handles);

% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in umaskLoad.
function umaskLoad_Callback(hObject, eventdata, handles)
% hObject    handle to umaskLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  c_dir= pwd;

   % if exist('img_history.mat','file')
   
   mname=mfilename('fullpath');
  dname=fileparts(mname);
  
   if exist(fullfile(dname,'umaskLoadExts.mat'),'file')
       load(fullfile(dname,'umaskLoadExts.mat'),'exts');
   else
       exts={'*.mat';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.hdr';'*.HEAD';'*.sdt'};
   end
   
  [files, dir_Name] = uigetfile(exts,'Load Underlay','MultiSelect','off');
   if ~dir_Name
        return;
    end
  [tmp,tmp2,ext]=fileparts(files);
  
  exts=unique(cat(1,{['*',ext]},exts),'stable');
  
  save(fullfile(dname,'umaskLoadExts.mat'),'exts');
  
  cd(c_dir);
  
    if ~iscell(files)
        temp = files;
        files = cell(1);
        files{1} = temp;
    end
    
    if length(files)>1
        error('select only one file');
    end
   
    
    setappdata(handles.figure1,'umfile',files{1});
    

    d=ri(fullfile(dir_Name,files{1}));
    d=double(d);
    setappdata(handles.figure1,'umdata',d);
    
    set(handles.showFull,'Value',true);
    
    dim=showImages(handles);
    
    set(handles.showFull,'Value',false);
    
    if length(dim)<4;
        dim(4)=1;
    end
  cbrik=get(handles.um_brik,'String');
  if str2num(cbrik)>dim(4)
      set(handles.um_brik,'String',num2str(dim(4)));
  end

  



function um_range_Callback(hObject, eventdata, handles)
% hObject    handle to um_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of um_range as text
%        str2double(get(hObject,'String')) returns contents of um_range as a double


% --- Executes during object creation, after setting all properties.
function um_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to um_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function um_brik_Callback(hObject, eventdata, handles)
% hObject    handle to um_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of um_brik as text
%        str2double(get(hObject,'String')) returns contents of um_brik as a double


% --- Executes during object creation, after setting all properties.
function um_brik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to um_brik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    rmappdata(handles.figure1,'umfile');
    rmappdata(handles.figure1,'umdata');
    
    showImages(handles);



function magn4roi_Callback(hObject, eventdata, handles)
% hObject    handle to magn4roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magn4roi as text
%        str2double(get(hObject,'String')) returns contents of magn4roi as a double
set(handles.showFull,'Value',true);
      
showImages(handles);

set(handles.showFull,'Value',false);


% --- Executes during object creation, after setting all properties.
function magn4roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magn4roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magn4roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in roiShiftAllSlices.
function roiShiftAllSlices_Callback(hObject, eventdata, handles)
% hObject    handle to roiShiftAllSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roiShiftAllSlices



function FlipROIVal_Callback(hObject, eventdata, handles)
% hObject    handle to FlipROIVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FlipROIVal as text
%        str2double(get(hObject,'String')) returns contents of FlipROIVal as a double


% --- Executes during object creation, after setting all properties.
function FlipROIVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FlipROIVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in logy.
function logy_Callback(hObject, eventdata, handles)
% hObject    handle to logy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logy


% --- Executes on button press in xdata.
function xdata_Callback(hObject, eventdata, handles)
% hObject    handle to xdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xdata



function rotationAngle_Callback(hObject, eventdata, handles)
% hObject    handle to rotationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotationAngle as text
%        str2double(get(hObject,'String')) returns contents of rotationAngle as a double
    udata=getappdata(handles.figure1,'udata');
    
    an=get(handles.rotationAngle,'String');
    an=str2num(an);
    for i=1:3
        if an(i)==0
            continue;
        end
        udata=rotate_images(udata,i,an(i));
    end
    
    setappdata(handles.figure1,'udata2',udata);
    
showImages(handles);

% --- Executes during object creation, after setting all properties.
function rotationAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveUnderlay.
function saveUnderlay_Callback(hObject, eventdata, handles)
% hObject    handle to saveUnderlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d=getappdata(handles.figure1,'udata');

[fname,dname]=uiputfile('*.mat','Underlay name');

voxsize=str2num(get(handles.u_voxsize,'String'));
center=str2num(get(handles.u_center,'String'));

if ~isempty(fname)
    save(fullfile(dname,fname),'d','center','voxsize');
end

function shiftVoxel_Callback(hObject, eventdata, handles)
% hObject    handle to shiftVoxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shiftVoxel as text
%        str2double(get(hObject,'String')) returns contents of shiftVoxel as a double


% --- Executes during object creation, after setting all properties.
function shiftVoxel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiftVoxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sliceStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiftVoxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliceStep_Callback(hObject, eventdata, handles)
% hObject    handle to text100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text100 as text
%        str2double(get(hObject,'String')) returns contents of text100 as a double

showImages(handles);

% --- Executes during object creation, after setting all properties.
function text100_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveOverLay.
function SaveOverLay_Callback(hObject, eventdata, handles)
% hObject    handle to SaveOverLay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

odata_rs=getappdata(handles.figure1,'odata_rs');

ind_u2o=getappdata(handles.figure1,'ind_u2o');
                
voxsize=str2num(get(handles.u_voxsize,'String'));
u_center=str2num(get(handles.u_center,'String'));
o_center=str2num(get(handles.o_center,'String'));


center=u_center;

nvox2shft = (u_center-o_center)./voxsize;

center(3)=(round(nvox2shft(3))-nvox2shft(3)).*voxsize(3)+o_center(3);


[fname,dname]=uiputfile('*.mat','Overlay name');

if ~isempty(fname)
    d=odata_rs;
    save(fullfile(dname,fname),'d','center','voxsize');
end

function odataShiftStep_Callback(hObject, eventdata, handles)
% hObject    handle to odataShiftStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of odataShiftStep as text
%        str2double(get(hObject,'String')) returns contents of odataShiftStep as a double


% --- Executes during object creation, after setting all properties.
function odataShiftStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to odataShiftStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function o_voxsize_Callback(hObject, eventdata, handles)
% hObject    handle to o_voxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of o_voxsize as text
%        str2double(get(hObject,'String')) returns contents of o_voxsize as a double


% --- Executes during object creation, after setting all properties.
function o_voxsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to o_voxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function o_center_Callback(hObject, eventdata, handles)
% hObject    handle to o_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of o_center as text
%        str2double(get(hObject,'String')) returns contents of o_center as a double


% --- Executes during object creation, after setting all properties.
function o_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to o_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function u_voxsize_Callback(hObject, eventdata, handles)
% hObject    handle to u_voxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of u_voxsize as text
%        str2double(get(hObject,'String')) returns contents of u_voxsize as a double


% --- Executes during object creation, after setting all properties.
function u_voxsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u_voxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function u_center_Callback(hObject, eventdata, handles)
% hObject    handle to u_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of u_center as text
%        str2double(get(hObject,'String')) returns contents of u_center as a double


% --- Executes during object creation, after setting all properties.
function u_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_resamp.
function calc_resamp_Callback(hObject, eventdata, handles)
% hObject    handle to calc_resamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(handles.figure1,'odata_rs',[]);
showImages(handles);


% --- Executes on button press in shift_overlay.
function shift_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to shift_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
step=get(handles.odataShiftStep,'String');
step=str2num(step);

odata_rs=getappdata(handles.figure1,'odata_rs');
odata_rs=circshift(odata_rs,step,3);
setappdata(handles.figure1,'odata_rs',odata_rs);

showImages(handles);


% --- Executes on button press in ROIMatchOverlay.
function ROIMatchOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to ROIMatchOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROIMatchOverlay
showImages(handles);

% --- Executes on button press in ClearROIPos.
function ClearROIPos_Callback(hObject, eventdata, handles)
% hObject    handle to ClearROIPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%setappdata(handles.figure1,'roi',[]);

udata=getappdata(handles.figure1,'udata');
odata_rs=getappdata(handles.figure1,'odata_rs');
ind_u2o = getappdata(handles.figure1,'ind_u2o');
slice=getappdata(gca,'slice');

%for i=min(1,1+ind_u2o):max(size(udata,3),size(odata_rs,3)-ind_u2o)
    i=slice;   
    j=1;
    while isappdata(handles.figure1,roi_pos_name(i,j))
      rmappdata(handles.figure1,roi_pos_name(i,j));
      j=j+1;
    end
    
    j=1;
    while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(i,j)))
      delete(findobj(handles.figure1,'Tag',roi_pos_name(i,j)));
      j=j+1;
    end
    
showImages(handles);



% --- Executes on button press in shiftRoiPosp.
function shiftRoiPosp_Callback(hObject, eventdata, handles)
% hObject    handle to shiftRoiPosp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ind_u2o=getappdata(handles.figure1,'ind_u2o');
ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');
if ~ROIMatchOverlay
    ind_u2o=0;
end


slice=getappdata(gca,'slice');

h=findobj(gca,'Tag',roi_pos_name(slice,1));

api = iptgetapi(h);

roi_pos = api.getPosition();
ud=get(handles.cropud,'String');
ud=str2num(ud);
lr=get(handles.croplr,'String');
lr=str2num(lr);

roi_pos(:,1)=roi_pos(:,1)+lr(1);

roi_pos(:,2)=roi_pos(:,2)+ud(1);

setappdata(gcf,roi_pos_name(slice+1,1),roi_pos);
delete(h);


showImages(handles);


    
    


% --- Executes on button press in shiftRoiPosm.
function shiftRoiPosm_Callback(hObject, eventdata, handles)
% hObject    handle to shiftRoiPosm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ind_u2o=getappdata(handles.figure1,'ind_u2o');
ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');
if ~ROIMatchOverlay
    ind_u2o=0;
end


slice=getappdata(gca,'slice');

h=findobj(gca,'Tag',roi_pos_name(slice,1));

api = iptgetapi(h);

roi_pos = api.getPosition();
ud=get(handles.cropud,'String');
ud=str2num(ud);
lr=get(handles.croplr,'String');
lr=str2num(lr);

roi_pos(:,1)=roi_pos(:,1)+lr(1);

roi_pos(:,2)=roi_pos(:,2)+ud(1);

setappdata(gcf,roi_pos_name(slice-1,1),roi_pos);
delete(h);


showImages(handles);


% --- Executes on button press in shiftVisibleROIs.
function shiftVisibleROIs_Callback(hObject, eventdata, handles)
% hObject    handle to shiftVisibleROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shiftVisibleROIs









% --- Executes on button press in detectuROI.
function detectuROI_Callback(hObject, eventdata, handles)
% hObject    handle to detectuROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ROIMatchOverlay,'Value',false);
calcROI_Callback(hObject, eventdata, handles);



% --- Executes on button press in PVSLabel.
function PVSLabel_Callback(hObject, eventdata, handles)
% hObject    handle to PVSLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROIMatchOverlay=get(handles.ROIMatchOverlay,'Value');

if ROIMatchOverlay
    roi=getappdata(gcf,'oroi');
else
    roi=getappdata(gcf,'uroi');
end

roi=merge_clusters_across_slices(roi);

if ROIMatchOverlay
    setappdata(gcf,'oroi',roi);
else
    setappdata(gcf,'uroi',roi);
end

showImages(handles);




    
    
    


% --- Executes on selection change in roiShapes.
function roiShapes_Callback(hObject, eventdata, handles)
% hObject    handle to roiShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiShapes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiShapes


% --- Executes during object creation, after setting all properties.
function roiShapes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
