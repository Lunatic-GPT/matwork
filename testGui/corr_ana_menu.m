function varargout = corr_ana_menu(varargin)
% CORR_ANA_MENU M-file for corr_ana_menu.fig
%      CORR_ANA_MENU, by itself, creates a new CORR_ANA_MENU or raises the existing
%      singleton*.
%
%      H = CORR_ANA_MENU returns the handle to a new CORR_ANA_MENU or the handle to
%      the existing singleton*.
%
%      CORR_ANA_MENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORR_ANA_MENU.M with the given input arguments.
%
%      CORR_ANA_MENU('Property','Value',...) creates a new CORR_ANA_MENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before corr_ana_menu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to corr_ana_menu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help corr_ana_menu

% Last Modified by GUIDE v2.5 21-May-2009 22:20:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @corr_ana_menu_OpeningFcn, ...
                   'gui_OutputFcn',  @corr_ana_menu_OutputFcn, ...
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


% --- Executes just before corr_ana_menu is made visible.
function corr_ana_menu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to corr_ana_menu (see VARARGIN)

% Choose default command line output for corr_ana_menu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes corr_ana_menu wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = corr_ana_menu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.data.sid=get(handles.sid,'String');


varargout{1} = handles.data;
delete(hObject);



function sid_Callback(hObject, eventdata, handles)
% hObject    handle to sid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sid as text
%        str2double(get(hObject,'String')) returns contents of sid as a double


% --- Executes during object creation, after setting all properties.
function sid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
%return;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drois_Callback(hObject, eventdata, handles)
% hObject    handle to drois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drois as text
%        str2double(get(hObject,'String')) returns contents of drois as a double


% --- Executes during object creation, after setting all properties.
function drois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


