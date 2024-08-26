function varargout = phase_correction_gui(varargin)
% PHASE_CORRECTION_GUI MATLAB code for phase_correction_gui.fig
%      PHASE_CORRECTION_GUI, by itself, creates a new PHASE_CORRECTION_GUI or raises the existing
%      singleton*.
%
%      H = PHASE_CORRECTION_GUI returns the handle to a new PHASE_CORRECTION_GUI or the handle to
%      the existing singleton*.
%
%      PHASE_CORRECTION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE_CORRECTION_GUI.M with the given input arguments.
%
%      PHASE_CORRECTION_GUI('Property','Value',...) creates a new PHASE_CORRECTION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phase_correction_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phase_correction_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phase_correction_gui

% Last Modified by GUIDE v2.5 29-Jan-2016 18:05:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phase_correction_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @phase_correction_gui_OutputFcn, ...
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

% --- Executes just before phase_correction_gui is made visible.
function phase_correction_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phase_correction_gui (see VARARGIN)

% Choose default command line output for phase_correction_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using phase_correction_gui.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes phase_correction_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phase_correction_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on slider movement.
function sliderPhase_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

value=get(hObject,'Value');

set(handles.phase,'String',num2str(value));

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function sliderPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Max',180);
set(hObject,'Min',-180);
set(hObject,'Value',0);
set(hObject,'SliderStep',[1/360,1/36]);


% --- Executes on slider movement.
function sliderTrnc_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTrnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

value=get(hObject,'Value');
trnc=get(handles.trnc,'String');
trnc=str2num(trnc);

trnc(1)=value;

set(handles.trnc,'String',num2str(trnc));

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function sliderTrnc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTrnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


set(hObject,'Max',1000);
set(hObject,'Min',0);
set(hObject,'Value',0);
set(hObject,'SliderStep',[0.001,0.01]);



function coord_Callback(hObject, eventdata, handles)
% hObject    handle to coord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord as text
%        str2double(get(hObject,'String')) returns contents of coord as a double

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function coord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OpenFile.
function OpenFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    c_dir= pwd;

   % if exist('img_history.mat','file')
   
  mname=mfilename('fullpath');
  dname=fileparts(mname);
  
   if exist(fullfile(dname,'imageSelectionExts.mat'),'file')
       load(fullfile(dname,'imageSelectionExts.mat'),'exts');
   else
       exts={'*.HEAD';'*.mat'};
   end
   
  [files, dir_Name] = uigetfile(exts,'MultiSelect','on');
   if ~dir_Name
        return;
    end
  [tmp,tmp2,ext]=fileparts(files);
  
  exts=unique(cat(1,{['*',ext]},exts),'stable');
  
  save(fullfile(dname,'imageSelectionExts.mat'),'exts');
  
  cd(c_dir);
  
  d=ri(fullfile(dir_Name,files));
  
  set(handles.coord,'String',num2str(floor((size(d(:,:,:,1)))/2)));
  
  setappdata(gcf,'data',d);
  set(handles.sliderTrnc,'Max',size(d,4));
  set(handles.sliderTrnc,'SliderStep',[1/size(d,4),10/size(d,4)]);
  
   set(gcf,'Name',files);
    
  setappdata(gcf,'newAxisScale',1);
  showSpectrum(handles);

  
  
  
  
  function showSpectrum(handles)
        
  
      newAxisScale=getappdata(gcf,'newAxisScale');
      if ~newAxisScale        
          xl=xlim;
          yl=ylim;
      end
      
  d=getappdata(gcf,'data');
  
  ijk=get(handles.coord,'String');
  ijk=str2num(ijk);
  
  ph=get(handles.phase,'String');
  ph=str2num(ph);
  
  trnc=get(handles.trnc,'String');
  trnc=str2num(trnc);
  
  baseline=get(handles.baseline,'String');
  baseline=str2num(baseline);
  
  d=squeeze(d(ijk(1),ijk(2),ijk(3),trnc(1)+1:end-trnc(2)));
  
d=d-mean(d(end-baseline+1:end));
fid=get(handles.showFid,'Value');



sw=str2num(get(handles.sw,'String'));

if ~fid 
    hold off;

fd=fft(d);
fd=fftshift(fd);


fd=fd*exp(1i*ph/180*pi);

fd=squeeze(fd);

centerF=str2num(get(handles.centerF,'String'));
f0=str2num(get(handles.f0,'String'));

x=centerF+(-length(fd)/2:length(fd)/2-1)*sw/length(fd)/f0;

plot(handles.axes1,x,flipud(real(fd)));   

xlim([min(x),max(x)]);
set(gca,'FontSize',10,'XDir','reverse');
xlabel('ppm');
else
    
x=(0:length(d)-1)/sw*1000;
plot(handles.axes1,x,real(d),'b');
hold on;
plot(handles.axes1,x,imag(d),'r');
xlim([0,max(x)]);
xlabel('ms');

end
    
 if ~newAxisScale   
    xlim(xl);
    ylim(yl);
 end

 
  setappdata(gcf,'newAxisScale',0);
  
function trnc_Callback(hObject, eventdata, handles)
% hObject    handle to trnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trnc as text
%        str2double(get(hObject,'String')) returns contents of trnc as a double
value=get(hObject,'String');
value=str2num(value);

set(handles.sliderTrnc,'Value',value(1));

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function trnc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phase_Callback(hObject, eventdata, handles)
% hObject    handle to phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phase as text
%        str2double(get(hObject,'String')) returns contents of phase as a double
value=get(hObject,'String');
set(handles.sliderPhase,'Value',str2num(value));

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baseline_Callback(hObject, eventdata, handles)
% hObject    handle to baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baseline as text
%        str2double(get(hObject,'String')) returns contents of baseline as a double

showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showFid.
function showFid_Callback(hObject, eventdata, handles)
% hObject    handle to showFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showFid

  setappdata(gcf,'newAxisScale',1);
showSpectrum(handles);

function centerF_Callback(hObject, eventdata, handles)
% hObject    handle to centerF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of centerF as text
%        str2double(get(hObject,'String')) returns contents of centerF as a double
showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function centerF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centerF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sw_Callback(hObject, eventdata, handles)
% hObject    handle to sw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sw as text
%        str2double(get(hObject,'String')) returns contents of sw as a double

showSpectrum(handles);
% --- Executes during object creation, after setting all properties.
function sw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0_Callback(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0 as text
%        str2double(get(hObject,'String')) returns contents of f0 as a double
showSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
