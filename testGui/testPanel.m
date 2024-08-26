function varargout = testPanel(varargin)
% TESTPANEL M-file for testPanel.fig
%      TESTPANEL, by itself, creates a new TESTPANEL or raises the existing
%      singleton*.
%
%      H = TESTPANEL returns the handle to a new TESTPANEL or the handle to
%      the existing singleton*.
%
%      TESTPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTPANEL.M with the given input arguments.
%
%      TESTPANEL('Property','Value',...) creates a new TESTPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testPanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testPanel

% Last Modified by GUIDE v2.5 07-Aug-2009 20:30:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @testPanel_OutputFcn, ...
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


% --- Executes just before testPanel is made visible.
function testPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testPanel (see VARARGIN)

% Choose default command line output for testPanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testPanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testPanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

value = get(handles.slider1,'Value');
pos=get(handles.uipanel1,'Position');
pos(2)=-0.5*value;
set(handles.uipanel1,'position',pos);
disp(get(handles.uipanel1,'Units'));
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


