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
% 10/5/2020: removed magn4roi box
%            added tooltips and removed labels for several boxes
% to do:
% Voxoffset, Voxstep, start trace: tooltips wrong


if ~verLessThan('MATLAB','9.7')
    set( groot , 'defaultAxesCreateFcn' , 'disableDefaultInteractivity(gca)' )
end
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
datacursormode off;
set(gcf,'Pointer','arrow');
% End initialization code - DO NOT EDIT

% --- Executes just before imageGallery is made visible.


function pos = myupdatefcn(obj,event_obj)
% Display 'Time' and 'Amplitude'
pos = get(event_obj,'Position');
plot_TimeCourse;


handles=getappdata(gcf,'handles');
ROIMatchOverlay = get(handles.ROIMatchOverlay,'Value');

if ROIMatchOverlay
    roi=getappdata(gcf,'oroi_rs');
    if isempty(roi)
       roi=getappdata(gcf,'oroi');
    end
else
    roi=getappdata(gcf,'uroi');
end

udata=getappdata(gcf,'udata');

odata=getappdata(gca,'odata_rs');

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
    str=sprintf('%s; o=%5.4f',str,odata(pos(1)+crop_lr(1),pos(2)+crop_ud(1)));
end

if ~isempty(roi)
    str=sprintf('%s; roi=%d',str,roi(pos(1)+crop_lr(1),pos(2)+crop_ud(1),slice,mbrik));
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
use_pwd=get(handles.pwd,'Value');

dataDir=get_dataDir(use_pwd,'underlay');


    h_udata=getappdata(handles.figure1,'udata');
    if isempty(h_udata)
        dataDir=pwd;
    end
    
if exist(fullfile(dname,'imageSelectionExts.mat'),'file')
    load(fullfile(dname,'imageSelectionExts.mat'),'exts');
else
    exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
end

[files, dataDir] = uigetfile(exts,'Load Underlay',dataDir,'MultiSelect','on');
if ~dataDir
    return;
end
[tmp,tmp2,ext]=fileparts(files);

exts=unique(cat(1,{['*',ext]},exts),'stable');

save(fullfile(dname,'imageSelectionExts.mat'),'exts');
save_dataDir(dataDir,'underlay');


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
setappdata(handles.figure1,'dname',dataDir);
%setappdata(handles.figure1,

if strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
    files_tmp=dir([dataDir,'/*',ext]);
    if length(files_tmp)>1
        ns=inputdlg(sprintf('Input Number of slices (1-%d)',length(files_tmp)));
    else
        ns={'1'};
    end
    d=ri(fullfile(dataDir,files{1}),str2num(ns{1}));
    d=double(d);
elseif strcmp(ext,'.mat')
    
    tmp=whos('-file',fullfile(dataDir,files{1}));
    
    for i=1:length(tmp)
        nm{i}=tmp(i).name;
    end
    
    if any(strcmp(nm,'voxsize'))
        voxsize= ri(fullfile(dataDir,files{1}),[],[],'voxsize');
    else
        voxsize=[1,1,1];
    end
    
    set(handles.u_voxsize,'String',sprintf('%3.1f,%3.1f,%3.1f',voxsize(:)));
    
    if any(strcmp(nm,'center'))
        center= ri(fullfile(dataDir,files{1}),[],[],'center');
    else
        center=[0,0,0];
    end
    
    set(handles.u_center,'String',sprintf('%3.1f,%3.1f,%3.1f',center));
    
    if any(strcmp(nm,'d'))
        d= ri(fullfile(dataDir,files{1}),[],[],'d');
    else
        
        d=ri(fullfile(dataDir,files{1}));
    end
    orient=get_orient(fullfile(dataDir,files{1}));
    setappdata(handles.figure1,'u_orient',orient);
    
elseif strcmp(ext,'.gz') 
    nii=load_untouch_niigz(fullfile(dataDir,files{1}));
    
    d=nii.img;
    nii=rmfield(nii,'img');

    setappdata(handles.figure1,'nii',nii);
    
elseif strcmp(ext,'.nii')   
    nii=load_untouch_nii(fullfile(dataDir,files{1}));
    d=nii.img;
     nii=rmfield(nii,'img');

    setappdata(handles.figure1,'nii',nii);
    
else
    
    d=ri(fullfile(dataDir,files{1}));
end

d=double(d);
setappdata(handles.figure1,'udata',d);
setappdata(handles.figure1,'odata_rs',[]); % overlay data may need to be resampled since underlay is changed.
setappdata(handles.figure1,'odata',[]);
setappdata(handles.figure1,'uroi',[]);
setappdata(handles.figure1,'oroi',[]);

setappdata(handles.figure1,'uroifile','');
setappdata(handles.figure1,'oroifile','');

if isappdata(handles.figure1,'udata2')
    rmappdata(handles.figure1,'udata2');
end

if isappdata(handles.figure1,'ImageTransforms')
    rmappdata(handles.figure1,'ImageTransforms');
end

set(handles.showFull,'Value',true);

dim=showImages(handles);

set(handles.showFull,'Value',false);

if length(dim)<4
    dim(4)=1;
end
set(handles.ts_brik,'String',sprintf('1 %d',dim(4)));
cbrik=get(handles.cur_brik,'String');
if str2num(cbrik)>dim(4)
    set(handles.cur_brik,'String',num2str(dim(4)));
end

if  strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
    [vox,center]=dcmDimCenter(dataDir);
    
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

function dataDir=get_dataDir(use_pwd,type)


if strcmp('underlay',type)
    fname='UnderlayDir.mat';
elseif strcmp('overlay',type)
    fname='OverlayDir.mat';
elseif strcmp('underlay mask',type)
     fname='UnderlayMaskDir.mat';
elseif strcmp('overlay mask',type)
     fname='OverlayMaskDir.mat';
else
    error('unknown type');
end
    

mname=mfilename('fullpath');
dname_img=fileparts(mname);

if exist(fullfile(dname_img,fname),'file') && ~use_pwd
    load(fullfile(dname_img,fname),'dataDir');
else
    dataDir=pwd;
end

function save_dataDir(d,type)

if strcmp('underlay',type)
    fname='UnderlayDir.mat';
elseif strcmp('overlay',type)
    fname='OverlayDir.mat';
elseif strcmp('underlay mask',type)
     fname='UnderlayMaskDir.mat';
elseif strcmp('overlay mask',type)
     fname='OverlayMaskDir.mat';
else
    error('unknown type');
end
    
dataDir=d;
mname=mfilename('fullpath');
dname_img=fileparts(mname);
if ~strcmp(d,pwd)
     save(fullfile(dname_img,fname),'dataDir');
end
    
% --- Executes on button press in load_roi.
function loadroi_Callback(hObject, eventdata, handles)
% hObject    handle to load_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mname=mfilename('fullpath');
dname_img=fileparts(mname);

use_pwd=get(handles.pwd,'Value');


roiMatchOverlay = get(handles.ROIMatchOverlay,'Value');

if ~roiMatchOverlay
dataDir = get_dataDir(use_pwd,'underlay mask');
else
dataDir = get_dataDir(use_pwd,'overlay mask');
    
end


if exist(fullfile(dname_img,'roiExts.mat'),'file')
    load(fullfile(dname_img,'roiExts.mat'),'exts');
else
    exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
end

[fname,dataDir]=uigetfile(exts,'Select the roi file',dataDir);

if ~dataDir
    return;
end
if ~roiMatchOverlay
    save_dataDir(dataDir,'underlay mask');
else
    save_dataDir(dataDir,'overlay mask');
end

[tmp,tmp2,ext]=fileparts(fname);

exts=unique(cat(1,{['*',ext]},exts),'stable');

save(fullfile(dname_img,'roiExts.mat'),'exts');



u_voxsize=str2num(get(handles.u_voxsize,'String'));
u_center=str2num(get(handles.u_center,'String'));

o_voxsize=str2num(get(handles.o_voxsize,'String'));
o_center=str2num(get(handles.o_center,'String'));

udata=getappdata(gcf,'udata');


if strcmp(ext,'.mat')
    
    tmp=load(fullfile(dataDir,fname));
    
    if (isfield(tmp,'center'))
        center=tmp.center;
    else
         if ~roiMatchOverlay
           center = u_center;
         else
           center = o_center;
         end
    end
    
    if (isfield(tmp,'voxsize'))
        voxsize=tmp.voxsize;
    else
         if ~roiMatchOverlay
          voxsize=u_voxsize;
         else
             voxsize=o_voxsize;
         end
    end
    

    roi=ri(fullfile(dataDir,fname));
    
    
else
    roi=ri(fullfile(dataDir,fname));
    
    if ~roiMatchOverlay
        center = u_center;
         voxsize=u_voxsize;
    else
        center = o_center;
         voxsize=o_voxsize;
    end
    
end

voxsize=voxsize(1:3);
u_voxsize=u_voxsize(1:3);

ignoreVoxSize =get(handles.ignoreVoxSize,'Value');


  roi=forwardImageTransform(roi,handles);
  
if roiMatchOverlay
    
    setappdata(handles.figure1,'oroi',roi);
    setappdata(handles.figure1,'oroifile',fullfile(dataDir,fname));
    setappdata(handles.figure1,'oroi_rs',[]);
else
    
    setappdata(handles.figure1,'uroifile',fullfile(dataDir,fname));
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

usepwd=get(handles.pwd,'Value');

dataDir=get_dataDir(usepwd,'overlay');

if exist(fullfile(dname_img,'overlayExts.mat'),'file')
    try
      load(fullfile(dname_img,'overlayExts.mat'),'exts');
    catch
      exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
    end
else
    exts={'*.hdr';'*.nii';'*.dcm';'*.IMA';'data-set_*.xml';'*.mat';'*.HEAD';'*.sdt'};
end

[files, dataDir]=uigetfile(exts,'Select the overlay file',dataDir);

if ~dataDir
    return;
end

[tmp,tmp2,ext]=fileparts(files);

exts=unique(cat(1,{['*',ext]},exts),'stable');

save(fullfile(dname_img,'overlayExts.mat'),'exts');

save_dataDir(dataDir,'overlay');

cd(c_dir);

if iscell(files)
    error('select only one file');
end

if ~dataDir
    return;
end

if  strcmp(ext,'.IMA') || strcmp(ext,'.dcm')
    [vox,center]=dcmDimCenter(dataDir);
    
    vox_str=sprintf('%5.3f,',vox);
    set(handles.o_voxsize,'String',vox_str(1:end-1));
    center_str=sprintf('%5.3f,',center);
    set(handles.o_center,'String',center_str(1:end-1));
    
end

setappdata(handles.figure1,'ofile',files);
setappdata(handles.figure1,'odir_Name',dataDir);

if strcmp(ext,'.mat')
    
    try
        voxsize= ri(fullfile(dataDir,files),[],[],'voxsize');
    catch
        voxsize=[1,1,1];
    end
    
        set(handles.o_voxsize,'String',sprintf('%3.1f,%3.1f,%3.1f',voxsize(:)));
        try
        center= ri(fullfile(dataDir,files),[],[],'center');
        catch
            center=[0,0,0];
        end
        set(handles.o_center,'String',sprintf('%3.1f,%3.1f,%3.1f',center(:)));
        
        d= ri_d1(fullfile(dataDir,files));
   try
        orient=get_orient(fullfile(dataDir,files));
   catch
       orient='';
   end
        setappdata(handles.figure1,'o_orient',orient);

    
else
    d=read_afni_sdt_images(fullfile(dataDir,files));
end

d=double(d);

d=forwardImageTransform(d,handles);
  
rng=min_max(d);
%   set(handles.edit12,'String',num2str(rng));
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


[files, dataDir] = uigetfile({'data-set_*.xml';'*.hdr';'*.mat';'*.HEAD';'*.sdt'},'MultiSelect','off');
cd(c_dir);

if ~dataDir
    return;
end

setappdata(handles.figure1,'mfile',files);
setappdata(handles.figure1,'mdir_Name',dataDir);
d=ri(fullfile(dataDir,files));

setappdata(handles.figure1,'mdata',d);
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

function checkbox2_Callback(hObject, eventdata, handles)
disp('');

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
    
    xr=xlim;
    yr=ylim;
    
    nd=get(handles.magn,'String');
    nd=str2double(nd);
    cm=getappdata(hAxes(ia),'colormap');
    slice=getappdata(hAxes(ia),'slice');
    if ~isempty(slice)&&any(slice==ssl)
        saveds(end+1)=slice;
    else
        continue;
    end
   
 %   ch=get(hAxes(ia),'Children');
    ch=findobj(hAxes(ia),'Type','image');
    c_image=get(ch,'CData');
  %  [c_image,cm]=add_roi_contour(c_image,cm,hAxes(ia),handles);

    xr(1)=ceil(xr(1));
    xr(2)=floor(xr(2));
    yr(1)=ceil(yr(1));
    yr(2)=floor(yr(2));
    
    if xr(1)<1
        xr(1)=1;
    end
    
    if yr(1)<1
        yr(1)=1;
    end
       
    if xr(2)>size(c_image,2)
        xr(2)=size(c_image,2);
    end
    
    if yr(2)>size(c_image,1)
        yr(2)=size(c_image,1);
    end
    
    
    c_image=c_image(yr(1):yr(2),xr(1):xr(2),:);
    
    y=repmat2(c_image,nd);
%     c_image(isnan(c_image))=0;
%     y=zeros(nd,size(c_image,1),nd,size(c_image,2),3);
%     for i=1:size(c_image,1)
%         for j=1:size(c_image,2)
%             ind=round(c_image(i,j));
%             if ind<1
%                 ind=1;
%             elseif ind>size(cm,1)
%                 ind=size(cm,1);
%             else
%             end
%             tmp=shiftdim(cm(ind,:),-3);
%             y(:,i,:,j,:)=repmat(tmp,[nd,1,nd,1,1]);
%         end
%     end
%     y=reshape(y,[size(c_image,1)*nd,size(c_image,2)*nd,3]);
    %  fname=sprintf('%s_%d_%d_%d_%dX%d.tif',strtok(get(handles.edit16,'String'),'.'),ud(1),ud(2),lr(1),lr(2),slice);
    fname=sprintf('%s_X%d.tif',strtok(get(handles.edit16,'String'),'.'),slice);
    imwrite(y,fname,'TIFF','Resolution',300);
    fprintf('%s saved\n',fname);
end
set(handles.saveslices,'String',num2str(saveds));
%set(hObject,'Enable','on');


function [img,cm]= add_roi_contour(img,cm,hAx,handles)


slice=getappdata(hAx,'slice');
img(img>size(cm,1))=size(cm,1);

roi=saveROI_Callback([],[], handles,false);

croplr=get(handles.croplr,'String');
cropud=get(handles.cropud,'String');
croplr=str2num(croplr);
cropud=str2num(cropud);

roi = roi(croplr(1)+1:end-croplr(2),cropud(1)+1:end-cropud(2),:);
roi=permute(roi,[2,1,3]);
if slice<=size(roi,3)
    gap_BW=bwmorph(roi(:,:,slice),'remove');
    
    img(gap_BW>0)=size(cm,1)+1;
    
    
    %cm=cat(1,cm,[1,1,1]);
    if sum(gap_BW(:))>0
        cm=cat(1,cm,[0,1,0]);
    end
    
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
%setappdata(handles.figure1,'uroifile',[]);

setappdata(handles.figure1,'oroi',[]);
setappdata(handles.figure1,'oroi_rs',[]);

udata=getappdata(handles.figure1,'udata');
ind_u2o=getappdata(handles.figure1,'ind_u2o');
odata=getappdata(handles.figure1,'odata');

smin=1+min([0,ind_u2o]);
smax=max([size(odata,3)-ind_u2o,size(udata,3)]);
noption=length(get(handles.roiShapes,'String'));

for i=smin:smax
    
    
    j=1;
    while isappdata(handles.figure1,roi_pos_name(i,j))
        rmappdata(handles.figure1,roi_pos_name(i,j));
        j=j+1;
    end
    
    for k=1:noption
        j=1;
        while ~isempty(findobj(handles.figure1,'Tag',roi_pos_name(i,j,k)))
            delete(findobj(handles.figure1,'Tag',roi_pos_name(i,j,k)));
            j=j+1;
        end
    end
    
    
end

showImages(handles);





function res=get_odata_rs_slice_center(handles)

o_center=str2num(get(handles.o_center,'String'));
u_voxsize=str2num(get(handles.u_voxsize,'String'));
u_center=str2num(get(handles.u_center,'String'));

center_u2o=u_center-o_center;

nvox2shft = (center_u2o./u_voxsize(1:3));


res=o_center(3)+(round(nvox2shft(3))-nvox2shft(3)).*u_voxsize(3);







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


% --- Executes on button press in detectROI.
function detectROI_Callback(hObject, eventdata, handles)
% hObject    handle to detectoROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


calcROI_Callback(hObject, eventdata, handles);


function calcROI_Callback(hObject, eventdata, handles)

sig=get(handles.sigma,'String');

sig=str2num(sig);

%%
matchOverlay=get(handles.ROIMatchOverlay,'Value');

if matchOverlay
    d=getappdata(handles.figure1,'odata');
    ind_u2o=getappdata(handles.figure1,'ind_u2o');
    brik=get(handles.obrik,'String');
    brik=str2num(brik);
     roi=getappdata(handles.figure1,'oroi');
    roifile=getappdata(handles.figure1,'oroifile');
else
    d=getappdata(handles.figure1,'udata');
    ind_u2o=0;
    brik=get(handles.cur_brik,'String');
    brik=str2num(brik);
     roi=getappdata(handles.figure1,'uroi');
    roifile=getappdata(handles.figure1,'uroifile');
end


if isempty(roi) || size(roi,1)~=size(d,1) || size(roi,2)~=size(d,2) || size(roi,3)~=size(d,3)
    roi=0*d(:,:,:,1);
end

sd_scale=get(handles.sigma,'String');
sd_scale=str2num(sd_scale);



vroi=get(handles.ROIVal,'String');
vroi=str2num(vroi);


    roi_new=saveROI_Callback(hObject, eventdata, handles,false);

    roi4det=roi_new>0&roi==0;
    
    [roi_ind,ind]=calc1ROI(d(:,:,:,brik),roi4det,sd_scale);


     roi2=roi*0;    
     roi2(ind)=roi_ind*vroi;     
     
thr=get(handles.ClusterThreshold);
thr=str2num(thr);

     roi2=clusterize2(roi2,thr,0);
     roi2 = roi2==1;
 
     
     vroi=unique(roi2(:));
     vroi(vroi==0)=[];
     

    
    dname=fileparts(roifile);
    save(fullfile(dname,'roi_before_calcROI'), 'roi');

    
    if matchOverlay
        setappdata(handles.figure1,'oroi',roi2);
    else
        setappdata(handles.figure1,'uroi',roi2);  
    end
showImages(handles);
    
%     set(handles.roiShiftAll,'Value',false);
%     roiShiftAll=false;
%     val=unique(vec(m(xl,yl,slice+ind_u2o)));
%     val(val==0)=[];
% 
% roival=get(handles.ROIVal,'String');
% roival=str2num(roival);
% 
% 
% %roival=ceil(roival/2)*2;
% %roic=(roi==roival-1);
% slice=getappdata(gca,'slice');
% 
% h=findobj(handles.figure1,'Tag',sprintf('roi_pos%d_1',slice));
% api = iptgetapi(h);
% roi_pos = api.getPosition();
% ud=get(handles.cropud,'String');
% ud=str2num(ud);
% lr=get(handles.croplr,'String');
% lr=str2num(lr);
% 
% roi_pos(:,1)=roi_pos(:,1)+lr(1);
% 
% roi_pos(:,2)=roi_pos(:,2)+ud(1);
% 
% roi_tmp=roipoly(0*d(:,:,1,1),roi_pos(:,2),roi_pos(:,1));
% 
% roi(:,:,slice+ind_u2o)=roi(:,:,slice+ind_u2o).*double(roi_tmp==0);
% 
% tmp=d(:,:,slice+ind_u2o,brik);
% 
% y=sort(tmp(roi_tmp>0));
% 
% if sig>0
%     mn=mean(y(round(end/8):end));
%     sd=std(y(round(end/8):end));
%     mask=(tmp>mn+sig*sd)&roi_tmp>0;
% else
%     mn=mean(y(1:round(end/8*7)));
%     sd=std(y(1:round(end/8*7)));
%     
%     mask=tmp<mn+sig*sd&roi_tmp>0;
% end
% %%
% mc=clusterize2(mask,2);
% roi(:,:,slice+ind_u2o)=roi(:,:,slice+ind_u2o)+mc.*mask;
% 
% if matchOverlay
%     setappdata(handles.figure1,'oroi',roi);
% else
%     setappdata(handles.figure1,'uroi',roi);
% end
% 
% 
function [roi_ind3,ind3]=calc1ROI(udata,roisd_scale)
% find voxesl above/below a threshold mean+/-sd_scal*sd,
% with mean and sd calculated from the ring outside of the ROI roi==val;
% output:
%     ind3 - the ROI containing after expanding the roi==val with a voxel
%     roi_ind3 - the mask values for the voxels defined by ind3.

ind=find(roi);

[ind2,ind3]=m_ext_calc(size(roi),ind);
tmp=roi(ind2);
ind2(tmp>0)=[];

mn=mean(udata(ind2));
sd=std(udata(ind2));

udata_tmp=udata(ind3);

if sd_scale>0
    roi_ind3=(udata_tmp>=mn+sd_scale*sd);
else
    roi_ind3=(udata_tmp<=mn+sd_scale*sd);
end

roi2=roi*0;
roi2(ind3)=roi_ind3;
roi2=clusterize2(roi2);
roi2=roi2==1;

roi_ind3=roi2(ind3);


function [ind2,ind3]=m_ext_calc(sz,ind)

% determine the outer boundary of an ROI which is defined by indices ind;
% output:
%   ind2: the outer boundary
%   ind3: the outer boundary and the inside
l=2;
ind2=zeros(length(ind)*(2*l+1)^3,1);
count=0;
for i=1:length(ind)
    
    ijk=ind2subb(sz,ind(i));
    
    if ijk(1)-l<1 || ijk(1)+l>sz(1) || ijk(2)-l<1 || ijk(2)+l>sz(2) ||ijk(3)-l<1 || ijk(3)+l>sz(3)
        continue;
    end
    
    for i1=-l:l
        for i2=-l:l
            for i3=-l:l            
                count=count+1;
                ind2(count)=ind(i)+i1+i2*sz(1)+i3*sz(1)*sz(2);
            end
        end
    end
    
end

ind3=ind2;
ind2=setdiff(ind2,ind);

ind2(ind2==0)=[];
ind3(ind3==0)=[];



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

val=uint16(val);
m=uint16(m);
if roiShiftAll
    
    m2=m(:,:,slice+ind_u2o);
else
    m2=0*m(:,:,slice+ind_u2o);
    
    for i=1:length(val)
        m2=m2+val(i)*uint16(m(:,:,slice+ind_u2o)==val(i));
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
elseif strcmp(tag,'roidown')
    m(:,:,slice+ind_u2o)=m3+circshift(m2,[0,step,0]);
else
    m(:,:,slice+ind_u2o)=m3+circshift(m2,[0,0,step]);
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
odata=getappdata(handles.figure1,'odata');

ind=get(handles.cur_brik,'String');
ind=str2num(ind);

oind=get(handles.obrik,'String');
oind=str2num(oind);

show_Overlay=get(handles.showOverlay,'Value');

if show_Overlay && ~isempty(odata)
    tmp=odata(:,:,:,oind);
else
    tmp=udata(:,:,:,ind);
end

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

it=getappdata(handles.figure1,'ImageTransforms');

it{end+1}=tag;

setappdata(handles.figure1,'ImageTransforms',it);

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


mdata=getappdata(handles.figure1,'mdata');
mdata=permute(mdata,new_order);
setappdata(handles.figure1,'mdata',mdata);

o_voxsize=str2num(get(handles.o_voxsize,'String'));
set(handles.o_voxsize,'String',num2str(o_voxsize(new_order(1:3))));

o_center=str2num(get(handles.o_center,'String'));
set(handles.o_center,'String',num2str(o_center(new_order(1:3))));

u_voxsize=str2num(get(handles.u_voxsize,'String'));
set(handles.u_voxsize,'String',num2str(u_voxsize(new_order(1:3))));

u_center=str2num(get(handles.u_center,'String'));
set(handles.u_center,'String',num2str(u_center(new_order(1:3))));


showImages(handles);




% --- Executes on button press in flipx_image.
function flip_image_Callback(hObject, eventdata, handles)
% hObject    handle to flipx_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tag=get(hObject,'Tag');
if strcmp(tag,'flipx_image')
    dim=1;
elseif strcmp(tag,'flipy_image')
    dim=2;
else
    dim=3;
end


it=getappdata(handles.figure1,'ImageTransforms');
it{end+1}=tag;
setappdata(handles.figure1,'ImageTransforms',it);

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


d=getappdata(handles.figure1,'mdata');
d=flip(d,dim);
setappdata(handles.figure1,'mdata',d);

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
val=get(handles.ZoomToVal,'String');
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
disp('');
showImages(handles);



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


function roi2=expand_roi(roi)

unq=unique(roi(:));
unq(unq==0)=[];
roi2=0*roi;
for i=1:length(unq)
    tmp=bwmorph(roi==unq(i),'dilate',5);
    roi2(tmp)=unq(i);
    
end




% --- Executes on button press in rmCon.
function rmCon_Callback(hObject, eventdata, handles)
% hObject    handle to rmCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



roi=getappdata(handles.figure1,'uroi');
sz=size(roi);



Draw1_Callback(hObject, eventdata, handles);



slice=getappdata(gca,'slice');
m=draw2roi(handles,sz,slice);
seed=zeros(sz);
seed(:,:,slice)=m;
roi=remove_con(roi,seed);
setappdata(handles.figure1,'uroi',roi);
showImages(handles);



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
    roi_name='oroi';
else
   roi_name='uroi';
end
 roi_orig=getappdata(handles.figure1,roi_name);

ROIVal=get(handles.ZoomToVal,'String');
ROIVal=str2num(ROIVal);

roi=roi_orig;

roi(roi~=ROIVal)=0;

%%
slice=getappdata(gca,'slice');
roi2(:,:,slice)=expand_roi(roi(:,:,slice));
setappdata(handles.figure1,roi_name,roi2);

showImages(handles);


pause(1);
%set(handles.showROI,'Value',cValue);

setappdata(handles.figure1,roi_name,roi_orig);

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

thr=get(handles.ClusterThreshold,'String');
thr=str2num(thr);

roi=clusterize2(roi>0,thr,0);

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
    
    
   % [c_image,cm]=add_roi_contour(c_image,cm,hAxes(ia),handles);
    
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



% --- Executes on button press in olleft.
function olshift_Callback(hObject, eventdata, handles)
% hObject    handle to olleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    step=get(handles.odataShiftStep,'String');
    step=str2num(step);
   % orient=getappdata(handles.figure1,'o_orient');
    
    odata=getappdata(handles.figure1,'odata');
   
    label=get(hObject,'String');
    
    
    if strcmp(label,'left')
        
        odata=circshift(odata,[-step,0,0,0]);
      %  orient.center(1)=orient.center(1)-orient.voxsize(1)*step;
        
        
    elseif strcmp(label,'right')
         odata=circshift(odata,[step,0,0,0]);
      %  orient.center(1)=orient.center(1)+orient.voxsize(1)*step;
    elseif strcmp(label,'up')
         odata=circshift(odata,[0,-step,0,0]);
      %  orient.center(2)=orient.center(2)-orient.voxsize(2)*step;
    elseif strcmp(label,'down')
      %  orient.center(2)=orient.center(2)+orient.voxsize(2)*step;
        
         odata=circshift(odata,[0,step,0,0]);
    elseif strcmp(label,'X')
         odata=circshift(odata,[0,0,step]);
      %  orient.center(3)=orient.center(3)+orient.voxsize(3)*step;
    end
    
   % set(handles.o_center,'String',sprintf('%3.1f,%3.1f,%3.1f',orient.center));
    
  % setappdata(handles.figure1,'o_orient',orient);
 setappdata(handles.figure1,'odata',odata);
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

[files, dataDir] = uigetfile(exts,'Load Underlay Mask','MultiSelect','off');
if ~dataDir
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


d=ri(fullfile(dataDir,files{1}));
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

if isappdata(handles.figure1,'umdata')
  rmappdata(handles.figure1,'umdata');
end

if isappdata(handles.figure1,'um')
rmappdata(handles.figure1,'umfile');
end
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

orient=getappdata(handles.figure1,'o_orient');


[fname,dname]=uiputfile('*.mat','Overlay name');
odata=getappdata(handles.figure1,'odata');
if ~isempty(fname)
    save(fullfile(dname,fname),'odata');
    
    if ~isempty(orient)
       mat_addv_instruct(fullfile(dname,fname),orient);
    end
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
setappdata(handles.figure1,'odata_rs',[]);
setappdata(handles.figure1,'oroi_rs',[]);

showImages(handles);

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
setappdata(handles.figure1,'odata_rs',[]);
setappdata(handles.figure1,'oroi_rs',[]);

showImages(handles);


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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listbox2.
function listbox2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ROISize_Callback(hObject, eventdata, handles)
% hObject    handle to text53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text53 as text
%        str2double(get(hObject,'String')) returns contents of text53 as a double


% --- Executes during object creation, after setting all properties.
function ROISize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ignoreVoxSize.
function ignoreVoxSize_Callback(hObject, eventdata, handles)
% hObject    handle to ignoreVoxSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ignoreVoxSize




% --- Executes during object creation, after setting all properties.
function cp_roi_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cp_roi_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cp_roi.
function cp_roi_Callback(hObject, eventdata, handles)
% hObject    handle to cp_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=getappdata(handles.figure1,'uroi');
sl_s=get(handles.cp_roi_source,'String');

sl_s=str2num(sl_s);

sl=getappdata(gca,'slice');


a(:,:,sl)=a(:,:,sl_s);

setappdata(handles.figure1,'uroi',a);
showImages(handles);



% --- Executes on button press in Connect.
function Connect_Callback(hObject, eventdata, handles)
% hObject    handle to Connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  roi=getappdata(handles.figure1,'uroi');
     
  val=get(handles.ROIVal,'String');
  val=str2double(val);
  
  a=clusterize2(roi==val);
  
  if max(a(:))>2
     disp('More than 2 clusters found; cannot proceed'); 
     return;
  end
  
  pos1=roiCOM(a==1);
  pos2=roiCOM(a==2);
  
  
   dist=sos(pos1-pos2,2);
   
                x=linspace(pos1(1),pos2(1),2*ceil(dist));
                y=linspace(pos1(2),pos2(2),2*ceil(dist));
                z=linspace(pos1(3),pos2(3),2*ceil(dist));
                
                tmp=0*roi;
                for ix=1:length(x)
                    tmp(round(x(ix)),round(y(ix)),round(z(ix)))=val;
                end
             roi(tmp>0)=tmp(tmp>0);
             

  setappdata(handles.figure1,'uroi',roi);


% --- Executes on button press in EraseVC.
function EraseVC_Callback(hObject, eventdata, handles)
% hObject    handle to EraseVC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uroi=getappdata(handles.figure1,'uroi');

Draw1_Callback(hObject, eventdata, handles);
[tmp,roi_draw]=saveROI_Callback(hObject, eventdata, handles,false);

tmp=unique(uroi(roi_draw>0));



for i=1:length(tmp)
    if tmp(i)==0
        continue;
    end
   uroi(uroi==tmp(i))=0;  
end

setappdata(handles.figure1,'uroi',uroi);
showImages(handles);



% --- Executes on selection change in listbox7.
function listbox7_Callback(hObject, eventdata, handles)
% hObject    handle to listbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox7


% --- Executes during object creation, after setting all properties.
function listbox7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function res=filename(a)

[res,name,suf]=fileparts(a);

if isempty(res)
    res=a;
else
    res=[name,suf];
end

function orient=get_orient(fname)


if strcmp(fname(end-2:end),'.gz') 
nii=load_untouch_niigz(fname);
orient=get_orient_from_nii(nii);

elseif strcmp(fname(end-3:end),'.nii')
nii=load_untouch_niigz(fname);
orient=get_orient_from_nii(nii);

else
orient=load(fname,'voxsize','center','rotmat','pos','orient');   
end





function ClusterThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ClusterThreshold as text
%        str2double(get(hObject,'String')) returns contents of ClusterThreshold as a double


% --- Executes during object creation, after setting all properties.
function ClusterThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClusterThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in addCluster.
function addCluster_Callback(hObject, eventdata, handles)
% hObject    handle to addCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% hObject    handle to rmCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



data=getappdata(handles.figure1,'udata');

sz=size(data);

roi=getappdata(handles.figure1,'uroi');
if isempty(roi)
    roi=zeros(sz);
end

global tmpROIVal;
global seedVal;
tmpROIVal=9999;
seedVal=10000;
seed=(seedVal==roi);
if ~any(seed(:))
    
    h=impoly(gca);
    set(h,'Tag','ClusterSeed');
    
    api = iptgetapi(h);
    roi_pos = api.getPosition();
    ud=get(handles.cropud,'String');
    ud=str2num(ud);
    lr=get(handles.croplr,'String');
    lr=str2num(lr);
    
    roi_pos(:,1)=roi_pos(:,1)+lr(1);
    
    roi_pos(:,2)=roi_pos(:,2)+ud(1);
    
    slice=getappdata(gca,'slice');
    seed=zeros(size(roi));
    seed(:,:,slice)=roipos2roi(roi_pos,sz(1:2));
    
end


roi(roi==tmpROIVal)=0;
thr=get(handles.Threshold,'String');
thr=str2double(thr);
roi(seed>0)=seedVal;

cluster_connected_with_seed(data>thr,seed,roi,handles);
%roi(clus>0)=tmpROIVal;

%setappdata(handles.figure1,'uroi',roi);
%showImages(handles);




function mask=remove_con(mask,seed)

tmp=cluster_connected_with_seed(mask,seed);
mask(tmp>0)=0;


function res=cluster_connected_with_seed(mask,seed,roi,handles)


res=mask;
a=find(seed>0);

cluster_temp=ind2subb(size(mask),a);

ivc=1;

nvc=length(a);
global tmpROIVal;
while ivc<=nvc
    
    nb = find_neighbors(cluster_temp(ivc,:),mask);
    
    nnb = size(nb,1);
    if nnb > 0
        cluster_temp(nvc+1:nvc+nnb,:) = nb;
        nvc=nvc+nnb;
        for a=1:nnb
            mask(nb(a,1),nb(a,2),nb(a,3))=0;
            if exist('roi','var')
              roi(nb(a,1),nb(a,2),nb(a,3))=tmpROIVal;
            end
        end
        if exist('roi','var')
           setappdata(handles.figure1,'uroi',roi);
           showImages(handles);
        end

    end
    ivc = ivc+1;
end

res(mask>0)=0;


function pwd_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double



function Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double


% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in acceptCluster.
function acceptCluster_Callback(hObject, eventdata, handles)
% hObject    handle to acceptCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



roi=getappdata(handles.figure1,'uroi');

roiVal=get(handles.figure1,'ROIVal');
roiVal=str2double(roiVal);
global tmpROIVal;
global seedVal;

roi(roi==tmpROIVal|roi==seedVal)=roiVal;



setappdata(handles.figure1,'uroi',roi);
showImages(handles);

  







% --- Executes during object creation, after setting all properties.
function thresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function thresholdRange_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdRange as text
%        str2double(get(hObject,'String')) returns contents of thresholdRange as a double


% --- Executes during object creation, after setting all properties.
function thresholdRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

   

% --- Executes on button press in traceFinish.
function traceFinish_Callback(hObject, eventdata, handles)
% hObject    handle to traceFinish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function searchRadius_Callback(hObject, eventdata, handles)
% hObject    handle to searchRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of searchRadius as text
%        str2double(get(hObject,'String')) returns contents of searchRadius as a double


% --- Executes during object creation, after setting all properties.
function searchRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to searchRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ZoomToVal_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomToVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZoomToVal as text
%        str2double(get(hObject,'String')) returns contents of ZoomToVal as a double


% --- Executes during object creation, after setting all properties.
function ZoomToVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZoomToVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveROIAs.
function saveROIAs_Callback(hObject, eventdata, handles)
% hObject    handle to saveROIAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Hypointense.
function Hypointense_Callback(hObject, eventdata, handles)
% hObject    handle to Hypointense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hypointense
