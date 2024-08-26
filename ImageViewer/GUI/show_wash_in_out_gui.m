function varargout = show_wash_in_out_gui(varargin)
% SHOW_WASH_IN_OUT_GUI M-file for show_wash_in_out_gui.fig
%      SHOW_WASH_IN_OUT_GUI, by itself, creates a new SHOW_WASH_IN_OUT_GUI or raises the existing
%      singleton*.
%
%      H = SHOW_WASH_IN_OUT_GUI returns the handle to a new SHOW_WASH_IN_OUT_GUI or the handle to
%      the existing singleton*.
%
%      SHOW_WASH_IN_OUT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_WASH_IN_OUT_GUI.M with the given input arguments.
%
%      SHOW_WASH_IN_OUT_GUI('Property','Value',...) creates a new SHOW_WASH_IN_OUT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_wash_in_out_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_wash_in_out_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_wash_in_out_gui

% Last Modified by GUIDE v2.5 07-Aug-2009 15:56:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_wash_in_out_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @show_wash_in_out_gui_OutputFcn, ...
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


% --- Executes just before show_wash_in_out_gui is made visible.
function show_wash_in_out_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_wash_in_out_gui (see VARARGIN)

% Choose default command line output for show_wash_in_out_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
setappdata(handles.figure1,'params',varargin{1});

params = varargin{1};
data_dir = get(params,'Select data directory'); 

temp=load(fullfile(data_dir,'detailed_analysis_results.mat'),'results_detailed_analysis');
res_d = temp.results_detailed_analysis;
temp=load(fullfile(data_dir,'analysis_results.mat'),'results','params');
fslice = get(temp.params,'First slice');
lslice = get(temp.params,'Total number of slices')+fslice -1;
res = temp.results;
vp = 0.625*0.625*2/1000;
v=res{end-4,1}*vp;
d = (6*v/pi)^(1/3); 
str1 = sprintf('Slice %d-%d:\n',fslice,lslice);
str2 = sprintf('%s %3.2f%% \n%s %3.2f cm3, eff. size %3.2f cm',res_d{end-2,2},100*res_d{end-2,1},res{end-4,2},v,d);
set(handles.text1,'String',[str1,str2],'FontSize',12);


% UIWAIT makes show_wash_in_out_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function show_images(handles,params)

nslices = get(params,'Total number of slices');
fslice =get(params,'First slice');
folder = get(params,'Select image directory');  % folder for downloaded images. Also the folder for saving the shifted images.
data_dir = get(params,'Select data directory'); 

setappdata(handles.figure1,'nslices',nslices);

h_unit = 1/min(6,nslices)*0.8;
w_unit = 0.17;
rPitch = 1/min(6,nslices);
setappdata(handles.figure1,'rPitch',rPitch);


%set(handles.uipanel1,'Units','Normalized','Position',[0,0,1,0.95]);
if nslices <=6  
    set(handles.slider1,'Enable','off');
else
    set(handles.slider1,'Enable','on','Value',1);
end


wi_arr=cell(1,nslices);
wo_arr=cell(1,nslices);
c_image_arr = cell(1,nslices);
roi_arr = cell(1,nslices);
for i=1:nslices
    
fname = sprintf('%02d_2.mat',i+fslice-1);
temp=load(fullfile(folder,fname),'image');

dir_str = dir(fullfile(data_dir,['*X',num2str(i+fslice-1)]));
temp2 = load(fullfile(data_dir,dir_str.name,'ROIs.mat'),'rect_pos','roi','tissue_BW_large','tissue2_BW_large','tissue2_BW_large_b');

rect_pos = temp2.tissue2_BW_large_b;
c_image = full2crop(temp.image,rect_pos);
roi = crop2full(temp2.roi,temp2.rect_pos,temp2.tissue_BW_large);
roi_arr{i} = full2crop(roi,rect_pos);

fname = sprintf('%02d_wash_in_norm_diff.mat',i+fslice-1);
wi = load(fullfile(folder, fname),'image');
wi = full2crop(wi.image,rect_pos);

fname = sprintf('%02d_wash_out.mat',i+fslice-1);
wo = load(fullfile(folder, fname),'image');
wo = full2crop(wo.image,rect_pos);

if ~exist('wi_min','var')
    wi_min = min(wi(:));
    wi_max = max(wi(:));
    img_max = max(c_image(:));
else
    if wi_min>min(wi(:))
      wi_min = min(wi(:));
    end
    if wi_max < max(wi(:))
        wi_max = max(wi(:));
    end
    if img_max < max(c_image(:))
        img_max = max(c_image(:));
    end
end

wi_arr{i} = wi;
wo_arr{i} = wo;
c_image_arr{i} = c_image; 

end

hAxes =zeros(nslices,4);
for i=1:nslices
% image1
hAxes(i,1) = axes('Parent',handles.uipanel1, 'position', [0, 1-rPitch*i, w_unit, h_unit]);
%hAxes = axes('Parent',handles.uipanel1);
 norm_image = 100*c_image_arr{i}/img_max;
image(ind2rgb2d(norm_image,gray(100)),'Parent',hAxes(i,1));
title(['Slice ',num2str(i+fslice-1),', t=2']);
axis off image;


% image 2
 norm_image = 100*c_image_arr{i}/img_max;
norm_image(norm_image<2) = 2;
%edge1 = bwmorph(roi,'remove') | bwmorph(imfill(tissue,'hole'),'remove') | bwmorph(imfill(tissue2,'hole'),'remove');
edge1 = bwmorph(roi_arr{i},'remove');

hAxes(i,2) = axes('Parent',handles.uipanel1, 'position', [0.25,1-rPitch*i,w_unit,h_unit]);
  
norm_image(edge1>0) =1;
cm = gray(100);
cm(1,:) = [0,1,0];
  
image(ind2rgb2d(norm_image,cm),'Parent',hAxes(i,2));
title(['Slice ',num2str(i+fslice-1),' with roi']);
axis off image;

% image 3
 
cm = jet(200);
cm(1,:) = 0;    

hAxes(i,3) = axes('Parent',handles.uipanel1,'position', [0.5,1-rPitch*i,w_unit,h_unit]);
  
wi_tmp = wi_arr{i}*100;
wi_tmp(wi_tmp<2) = 2;
wi_tmp(edge1>0) = 1;
rgb = ind2rgb2d(wi_tmp,cm);

    image(rgb,'Parent',hAxes(i,3));
    title(['Slice ',num2str(i+fslice-1),': Wash-in (percent change)']);
    
    
    
    colormap(jet(90));
    colorbar('YTickLabel',{'0','50','100','150','200'},'YTick',1:89/4:90);
    axis off image;
   set(hAxes(i,3), 'position', [0.5,1-rPitch*i,w_unit,h_unit]);
% image 4    
    jet_inv = flipud(jet(180));
    jet_inv(1,:)=0;
    hAxes(i,4) = axes('Parent',handles.uipanel1,'position', [0.75,1-rPitch*i,w_unit,h_unit]);
    
    wo_tmp = wo_arr{i}+91;
    wo_tmp(wo_tmp<2)=2;
    wo_tmp(edge1>0)=1;
    rgb = ind2rgb2d(wo_tmp,jet_inv);
    image(rgb,'Parent',hAxes(i,4));
    title(['Slice ',num2str(i+fslice-1),': Wash-out']);
    colorbar('YTickLabel',{'90','45','0','-45','-90'},'YTick',1:89/4:90);
    axis off image;
    set(hAxes(i,4), 'position', [0.75,1-rPitch*i,w_unit,h_unit]);
    
end

setappdata(handles.figure1,'hAxes',hAxes);

% --- Outputs from this function are returned to the command line.
function varargout = show_wash_in_out_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

params = getappdata(handles.figure1,'params');



show_images(handles,params);


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = get(hObject,'Value');
nslices = getappdata(handles.figure1,'nslices');
hAxes = getappdata(handles.figure1,'hAxes');
rPitch = getappdata(handles.figure1,'rPitch');
top = rPitch*nslices+val*(1-rPitch*nslices);

   for i=1:4
       for j=1:nslices
           pos = get(hAxes(j,i),'position');
            
           
           pos(2) = top-rPitch*j;
           
           set(hAxes(j,i),'position',pos);
       end
   end
           
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


