function varargout = DICOM_view(varargin)

% Copyright © 2007 Michigan State University
%    
% All rights are reserved.  Reproduction or transmission in whole or in part, in any form or
% by any means, electronic, mechanical or otherwise, is prohibited without the prior written
% consent of the copyright owner.
% 
% Filename:  DICOM_view
% Created:  12/29/2006
% Author:  Tobias Hahn
% Uses some constants ('rows', 'rPitch', 'cPitch' etc.) from an open
% .m-file created by Yonathan Nativ, Ben Gurion University
% 
% Description:  Displays a series of DICOM images; THIS IS THE NEW VERSION
% 
% Usage:
%  DICOM_view
%
% DICOM_VIEW M-file for DICOM_view.fig
%      DICOM_VIEW, by itself, creates a new DICOM_VIEW or raises the existing
%      singleton*.
%
%      H = DICOM_VIEW returns the handle to a new DICOM_VIEW or the handle to
%      the existing singleton*.
%
%      DICOM_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICOM_VIEW.M with the given input arguments.
%
%      DICOM_VIEW('Property','Value',...) creates a new DICOM_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DICOM_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DICOM_view_OpeningFcn via varargin.
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DICOM_view_OpeningFcn, ...
                   'gui_OutputFcn',  @DICOM_view_OutputFcn, ...
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


% --- Executes just before DICOM_view is made visible.
function DICOM_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DICOM_view (see VARARGIN)


% Choose default command line output for DICOM_view
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes DICOM_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DICOM_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global dirName 

% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)

P = fileparts(mfilename('fullpath'));
global dirName
dirName = uigetdir(P, 'Select DICOM directory (This is the folder where you find the different DICOM-series, sorted in folders with a max. of 25 images)');
if ~ischar(dirName)
    disp('no valid Directory selected.')
    return;
end

cd(dirName);

handles.dfolder=dirName;
guidata(hObject, handles);
SetFolder(handles, dirName);

% hObject    handle to imageSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetFolder (handles, dirName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prompt to determine the designated image

user_entry = input('Do you want to use an existing image list (type 1) or create a new one (type 2)?');
if user_entry == 1
    
    %% THIS WOULD BE WONDERFUL TO INCLUDE, BUT MATLAB HAS SOME ERROR HERE
    %load filename_array
    %%
     prompt = {'image#'};
     defans={'1'};
     fields = {'number'};
%     info=inputdlg(prompt, 'Designate image number:',1,defans);
%     if ~isempty(info) % see if user hit 'cancel'
%         info = cell2struct(info, fields);
%         myimage = str2num(info.number);
%     end

    info=inputdlg(prompt, 'Input the looked-up corresponding Acqu. Number',1,defans);
    if ~isempty(info) % see if user hit 'cancel'
        info = cell2struct(info, fields);
        i=1;
        filename_array(i,2) = str2num(info.number);
    end
    
    
    %image_found = 0; % flag to show if image is found
    %for i=1:size(filename_array,1)
        %if filename_array(i,1) == myimage
            image_found = 1;
            imagenumber = 1;
            index = 1;
            cols = 4;
            rows = ceil(imagenumber/cols);
            ratio = rows/cols;
            hight = ratio*0.95;
            pos = [0 -(hight-0.95) 0.95 hight];
            set(handles.uipanel1,'position',pos);
            set(handles.slider1,'enable','on','value',1);

            rPitch = 0.98/rows;
            cPitch = 0.98/cols;
            axWidth = 0.9/(cols);
            axHight = 0.9/(rows);

            hAxes = getappdata(handles.figure1,'hAxes');
            if ~isempty(hAxes)
                f = find ( ishandle(hAxes) & hAxes);
                delete(hAxes(f));
            end

            hAxes = zeros(imagenumber,1);
            [r c] = ind2sub([rows cols],index);
            x = 0.98-(c)*cPitch;
            y = 0.98-(r)*rPitch;
            hAxes(index) = axes( 'position', [x y axWidth axHight], 'Parent', handles.uipanel1);
            iname = num2str(filename_array(i,2));
            img = dicomread([handles.dfolder '\' iname]);
            plotImInAxis(img,hAxes(index));
        %end
    %end
    if image_found ==1
        setappdata(handles.figure1,'hAxes',hAxes);
    else
        disp('The image could not be found');
    end
elseif user_entry == 2
     if exist('filename_array.mat','file') ~= 0
         user_entry = input('There is already a filename array. Do you want to delete it (yes = 1) or move it to the MATLAB folder? (yes = 2)?');
         if user_entry == 1
            delete('filename_array.mat');
         elseif user_entry == 2
             movefile('filename_array.mat','C:\Documents and Settings\hahntobi\My Documents\MATLAB\filename_arrays')
         end
     end
    dfiles=dir(handles.dfolder);
    dfiles=dfiles(3:end);                   % avoid . and ..
    nfiles=length(dfiles); % number of images
    if nfiles<1
        disp('No DICOM images in this folder.')
        return;
    end
    filename_array = zeros(nfiles,2);
    index = 1;
    for i = 1:nfiles
        iname = dfiles(i).name; %the image name
        if isdicom([handles.dfolder '\' iname])
            info = dicominfo([handles.dfolder '\' iname]);
            imac = info.InstanceNumber;
            filename_array(index,1)=imac;
            filename_array(index,2)=str2num(iname); % WORKS ONLY IN THIS CASE SINCE FILES ARE NUMBERED       
            index = index + 1;
        end
    end
    %cd('C:\Documents and Settings\hahntobi\My Documents\MATLAB\filename_arrays\current filename_array');
    save 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\filename_arrays\current filename_array\filename_array' filename_array;
elseif user_entry == 3
    %%% THIS IS AT THE MOMENT ONLY FOR TEMPORARY USE - THEREFORE THE IMAGES
    %%% ARE PREDEFINED
    
    filename_array(1,2) = 996527;

    filename_array(2,2) = 99654;

    image_found = 1;
    imagenumber = 2;
    index = 1;
    cols = 4;
    rows = ceil(imagenumber/cols);
    ratio = rows/cols;
    hight = ratio*0.95;
    pos = [0 -(hight-0.95) 0.95 hight];
    set(handles.uipanel1,'position',pos);
    set(handles.slider1,'enable','on','value',1);

    rPitch = 0.98/rows;
    cPitch = 0.98/cols;
    axWidth = 0.9/(cols);
    axHight = 0.9/(rows);

    hAxes = getappdata(handles.figure1,'hAxes');
    if ~isempty(hAxes)
        f = find ( ishandle(hAxes) & hAxes);
        delete(hAxes(f));
    end

    hAxes = zeros(imagenumber,1);
    for i = 1:2
    [r c] = ind2sub([rows cols],i);
    x = 0.98-(c)*cPitch;
    y = 0.98-(r)*rPitch;
    hAxes(i) = axes( 'position', [x y axWidth axHight], 'Parent', handles.uipanel1);
    iname = num2str(filename_array(i,2));
    img = dicomread([handles.dfolder '\' iname]);
    plotImInAxis(img,hAxes(i));
    end
    
    if image_found ==1
        setappdata(handles.figure1,'hAxes',hAxes);
    else
        disp('The image could not be found');
    end
        
end

%% VERSION 2
% prompt = {'image#'};
% defans={'1'};
% fields = {'number'};
% info=inputdlg(prompt, 'Designate image number:',1,defans);
% if ~isempty(info) % see if user hit 'cancel'
%     info = cell2struct(info, fields);
%     myimage = str2num(info.number);
% end
% 
% dfiles=dir(handles.dfolder);
% dfiles=dfiles(3:end);                   % avoid . and ..
% nfiles=length(dfiles); % number of images
% if nfiles<1
%     disp('No DICOM images in this folder.')
%     return;
% end

% flag = 0;
% index = 1;
% imagenumber = 5;
% for i=1:nfiles 
%     iname = dfiles(i).name; %the image name
%     info = dicominfo([handles.dfolder '\' iname]);
%     imac = info.InstanceNumber;
%     if imac == myimage || imac == myimage + 224 || imac == myimage + 2*224 || imac == myimage + 3* 224|| imac == myimage + 4*224|| imac == myimage + 5*224
%         cols = 4;
%         rows = ceil(imagenumber/cols);
%         ratio = rows/cols;
%         hight = ratio*0.95;
%         pos = [0 -(hight-0.95) 0.95 hight];
%         set(handles.uipanel1,'position',pos);
%         set( handles.slider1,'enable','on','value',1);
% 
%         rPitch = 0.98/rows;
%         cPitch = 0.98/cols;
%         axWidth = 0.9/(cols);
%         axHight = 0.9/(rows);
% 
%         hAxes = getappdata(handles.figure1,'hAxes');
%         if ~isempty(hAxes)
%             f = find ( ishandle(hAxes) & hAxes);
%             delete(hAxes(f));
%         end
%         
%        
%         hAxes = zeros(imagenumber,1);
%         [r c] = ind2sub([rows cols],index);
%         x = 0.98-(c)*cPitch;
%         y = 0.98-(r)*rPitch;
%         hAxes(index) = axes( 'position', [x y axWidth axHight], 'Parent', handles.uipanel1);
%         img      = dicomread([handles.dfolder '\' iname]);
%         plotImInAxis(img,hAxes(index));
%         index = index + 1;
%         if index == (imagenumber+1)
%             flag = 1;
%         end
%     end
%     if flag ==1 
%         return;
%     end
% end
%% END VERSION 2

%% VERSION 1
% %Prompt for input of the designated series
% prompt={'series'};
% defans={'1'};
% fields = {'series'};
% info = inputdlg(prompt, 'Designate series:', 1, defans);
% if ~isempty(info)              %see if user hit cancel
%    info = cell2struct(info,fields);
%    myseries = str2num(info.series);  %convert string to number
% end
% 
% %fprintf('There are %d series available. ',c.seriesCount);
% 
% user_entry_series = myseries;
% 
% %user_entry_series = input('Please specify the designated series: ');
% 
% dfiles=dir(handles.dfolder);
% dfiles=dfiles(3:end);                   % avoid . and ..
% 
% % nfiles=length(dfiles); % number of images
% % if nfiles<1
% %     disp('No DICOM images in this folder.')
% %     return;
% % end
% 
% flag = 0;
% flag1 = 0;
% imagecounter = 0;
% 
% for i=1:length(dfiles) %check all folders for the first image of the wanted series
%     if flag1 == 1 
%             break;
%     end
%     foldername = dfiles(i).name;
%     dicomfiles = dir([handles.dfolder '\' foldername]);
%     dicomfiles = dicomfiles(3:end);
%     for j=1:length(dicomfiles) %check all dicomfiles in each folder
%         if flag1 == 1 % if we have already reached the last image of the series of interest
%             break;
%         end
%         fname = dicomfiles(j).name;
%         info = dicominfo([handles.dfolder '\' foldername '\' fname]);
%         series = info.SeriesNumber;
%         imac = info.ImagesInAcquisition;
%         if series == user_entry_series % if this is the correct series
%             flag = 1;
%             imagecounter = imagecounter + 1;
%             if imagecounter == 1 % if it's the first image of the series
%                 nfiles=imac;
%                 cols = 4;
%                 rows = ceil(nfiles/cols);
%                 ratio = rows/cols;
%                 hight = ratio*0.95;
%                 pos = [0 -(hight-0.95) 0.95 hight];
%                 set(handles.uipanel1,'position',pos);
%                 set( handles.slider1,'enable','on','value',1);
% 
%                 rPitch = 0.98/rows;
%                 cPitch = 0.98/cols;
%                 axWidth = 0.9/(cols);
%                 axHight = 0.9/(rows);
%                 
%                 hAxes = getappdata(handles.figure1,'hAxes');
%                 if ~isempty(hAxes)
%                     f = find ( ishandle(hAxes) & hAxes);
%                     delete(hAxes(f));
%                 end
%     
% 
%                 hAxes = zeros(nfiles,1);
%             end
%             [r c] = ind2sub([rows cols],imagecounter);
%             x = 0.98-(c)*cPitch;
%             y = 0.98-(r)*rPitch;
%             hAxes(imagecounter) = axes( 'position', [x y axWidth axHight], 'Parent', handles.uipanel1);
%             img      = dicomread([handles.dfolder '\' foldername '\' fname]);
%             plotImInAxis(img,hAxes(imagecounter));
%         elseif flag == 1
%                 flag1 = 1;
%         end
%     end
% end
%% END VERSION 1

%% IMPORTANT IF WE WANT TO USE VERIONS 1 or 2
% setappdata(handles.figure1,'hAxes',hAxes);
%%

% for i=1:nfiles
%     fname = dfiles(i).name;
%     img      = dicomread([handles.dfolder '\' fname]);
%     plotImInAxis(img);
% end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
    pos = get(handles.uipanel1,'position');
    hight = pos(4);
    if hight > 1
        val = get(hObject,'value');
        yPos = -val * (hight-0.95);
        pos(2) = yPos;
        set(handles.uipanel1,'position',pos);
    end
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function plotImInAxis(im, hAx)
    imageProp = { 'ButtonDownFcn'};
    imageVal = {'openSpecificImage'};
    imagesc(im, imageProp,imageVal,'parent',hAx );
    colormap(gray);
    axis(hAx,'image');
    axis(hAx,'off');    
    drawnow;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

P = fileparts(mfilename('fullpath'));
global dirName
dirName = uigetdir(P, 'Select DICOM directory (This is the folder where you find the different DICOM-series, sorted in folders with a max. of 25 images)');
if ~ischar(dirName)
    disp('no valid Directory selected.')
    return;
end

cd(dirName);

% handles.dfolder=dirName;
% guidata(hObject, handles);
% SetFolder(handles, dirName);


disp('Checking the folders for DICOM files. This may take time. Please wait.');
pause(1);

a = StudyCollection(); %Creates an empty StudyCollection object

a = Merge(a, dirName);

b = Inspect(a);

c = Inspect(b.study{1});

str1 = 'There are ';
str2 = num2str(c.seriesCount);
str3 = ' series available.';

string=[str1 str2 str3];

disp(string);




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


