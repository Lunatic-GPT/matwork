function varargout = quan(varargin)
% QUAN M-file for quan.fig
%      QUAN, by itself, creates a new QUAN or raises the existing
%      singleton*.
%
%      H = QUAN returns the handle to a new QUAN or the handle to
%      the existing singleton*.
%
%      QUAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUAN.M with the given input arguments.
%
%      QUAN('Property','Value',...) creates a new QUAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quan

% Last Modified by GUIDE v2.5 01-Mar-2008 11:55:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quan_OpeningFcn, ...
                   'gui_OutputFcn',  @quan_OutputFcn, ...
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


% --- Executes just before quan is made visible.
function quan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quan (see VARARGIN)

% Choose default command line output for quan
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = quan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function editText_Quan_Callback(hObject, eventdata, handles)


function editText_Quan_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton1_Callback(hObject, eventdata, handles)




