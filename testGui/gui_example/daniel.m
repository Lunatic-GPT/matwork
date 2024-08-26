function varargout = daniel(varargin)
% DANIEL M-file for daniel.fig
%      DANIEL, by itself, creates a new DANIEL or raises the existing
%      singleton*.
%
%      H = DANIEL returns the handle to a new DANIEL or the handle to
%      the existing singleton*.
%
%      DANIEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DANIEL.M with the given input arguments.
%
%      DANIEL('Property','Value',...) creates a new DANIEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before daniel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to daniel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help daniel

% Last Modified by GUIDE v2.5 01-Mar-2008 11:56:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @daniel_OpeningFcn, ...
                   'gui_OutputFcn',  @daniel_OutputFcn, ...
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


% --- Executes just before daniel is made visible.
function daniel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to daniel (see VARARGIN)

% Choose default command line output for daniel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes daniel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = daniel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function editText_Daniel_Callback(hObject, eventdata, handles)


function editText_Daniel_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function pushbutton1_Callback(hObject, eventdata, handles)
