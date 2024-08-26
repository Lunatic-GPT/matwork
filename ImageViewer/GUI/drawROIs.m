function varargout = drawROIs(varargin)
% DRAWROIS M-file for drawROIs.fig
%      DRAWROIS, by itself, creates a new DRAWROIS or raises the existing
%      singleton*.
%
%      H = DRAWROIS returns the handle to a new DRAWROIS or the handle to
%      the existing singleton*.
% 
% 
% This GUI allows you to open several images and batch process all of them
% 

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawROIs_OpeningFcn, ...
                   'gui_OutputFcn',  @drawROIs_OutputFcn, ...
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

% --- Executes just before drawROIs is made visible.
function drawROIs_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for drawROIs

handles.output = hObject;
%setappdata(handles.figure1,'actionItems',[handles.imageSelection]);

% Update handles structure
guidata(hObject, handles);

setappdata(handles.figure1,'o_image',varargin{1});
setappdata(handles.figure1,'colormap',varargin{2});

setappdata(handles.figure1,'rect_pos',varargin{3});
setappdata(handles.figure1,'o_hAx',varargin{4});
%setappdata(handles.figure1,'scale',varargin{5});

[dname,fname]=fileparts(varargin{5});
setappdata(handles.figure1,'root',fileparts(dname));
setappdata(handles.figure1,'fname',fname);


showImage(handles);

o_hAx = getappdata(handles.figure1,'o_hAx');
h=findobj(o_hAx,'Tag','inner_pos');
if isempty(h)
%    set(handles.save,'Enable','off');
    set(handles.saveas,'Enable','off');
    set(handles.clear,'Enable','off');
%    set(handles.autoroi,'Enable','off');
%    set(handles.autoroi2,'Enable','off');
end



% UIWAIT makes drawROIs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function showImage(handles)
    
    
    set(handles.figure1,'Units','pixels');
    set(handles.figure1,'Position',[913,307,653,574]);
    o_image=getappdata(handles.figure1,'o_image');
    rect_pos = getappdata(handles.figure1,'rect_pos');

    o_hAx = getappdata(handles.figure1,'o_hAx');

    axesProp = {'dataaspectratio' ,...
                                'Parent',...
                                'PlotBoxAspectRatio', ...
                                'xgrid' ,...
                                'ygrid'};
    axesVal = {[1,1,1] , ...
                            handles.uipanel1,...
                            [1 1 1]...
                            'off',...
                            'off'};
                        
    hAxes = axes( 'position', [0 0 1 0.9],axesProp,axesVal);
    setappdata(handles.figure1,'axes',hAxes);
    cm=getappdata(handles.figure1,'colormap');
    imshow(o_image,cm,'Parent',hAxes);
    
    
    axis(hAxes,'image');
    axis(hAxes,'off');
    slice=getappdata(o_hAx,'slice');
    title(sprintf('Slice %d',slice));
 %   fname2title(fname);
    
    h_par = get(get(o_hAx,'Parent'),'Parent');
    h_roi = findobj(o_hAx,'Tag','roi_pos');
   roi = getappdata(o_hAx,'roi');     
       
        
   if ~isempty(h_roi)
     api=iptgetapi(h_roi);
     roi_pos=api.getPosition();
  
  %   roi_pos(:,1) = roi_pos(:,1)-rect_pos(1)+1;
  %   roi_pos(:,2) = roi_pos(:,2)-rect_pos(2)+1;
  
      h=impoly(hAxes,roi_pos);
      set(h,'Tag','roi_pos');
      api = iptgetapi(h);
      api.setColor('green');
  
  
   end


  h = findobj(o_hAx,'Tag','inner_pos');
  if ~isempty(h)
   api=iptgetapi(h);
   inner_pos=api.getPosition();
   h = findobj(o_hAx,'Tag','outer_pos');
   api=iptgetapi(h);
   outer_pos=api.getPosition();
      
   inner_pos(:,1) = inner_pos(:,1)-rect_pos(1)+1;
   inner_pos(:,2) = inner_pos(:,2)-rect_pos(2)+1;
   h=impoly(hAxes,inner_pos);
   set(h,'Tag','inner_pos');
   api = iptgetapi(h);
   api.setColor('red');

   outer_pos(:,1) = outer_pos(:,1)-rect_pos(1)+1;
   outer_pos(:,2) = outer_pos(:,2)-rect_pos(2)+1;
   h=impoly(hAxes,outer_pos);
   set(h,'Tag','outer_pos');
   api = iptgetapi(h);
      api.setColor('red');
  end 

% --- Outputs from this function are returned to the command line.
function varargout = drawROIs_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
  
% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
%     imtool close all
    
  
           
    close(handles.figure1);
    
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in draw.
function draw_Callback(hObject, eventdata, handles)
% hObject    handle to draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hAx = getappdata(handles.figure1,'axes');


clear_Callback(hObject, eventdata, handles);

%title('Please select 1) INNER BOUND 2) OUTER BOUND 3) rough ROI');
 title('Please select 1) INNER BOUND 2) OUTER BOUND');
 
                inner = impoly(hAx,[]); %create mask
                api_inner = iptgetapi(inner);
                api_inner.setColor('red');
                set(inner,'Tag','inner_pos');
                
                outer = impoly(hAx,[]); %create mask
                api = iptgetapi(outer);
                api.setColor('red');
                set(outer,'Tag','outer_pos');
                
              
                               
                h = impoly(hAx,api_inner.getPosition);
                 api = iptgetapi(h);
                api.setColor('green');
             
                set(h,'Tag','roi_pos');
                
set(handles.clear,'Enable','on'); 
set(handles.autoroi,'Enable','on');
set(handles.autoroi2,'Enable','on');
%set(handles.save,'Enable','on');
set(handles.saveas,'Enable','on');

%fname2title(getappdata(handles.figure1,'fname'));


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hAxes = getappdata(handles.figure1,'axes');
%{
c_image = getappdata(handles.figure1,'c_image');

scale = getappdata(handles.figure1,'scale');
ch = findobj(hAxes,'Type','image');
colormap(gray(200));
set(ch,'CData',c_image);
caxis(scale);
%}
    c_image=getappdata(handles.figure1,'o_image');

%    setappdata(handles.figure1,'c_image',c_image);
    
    cm=getappdata(handles.figure1,'colormap');
    imshow(c_image,cm,'Parent',hAxes);
    
    
    axis(hAxes,'image');
    axis(hAxes,'off');
 %   fname2title(fname);


if (isappdata(handles.figure1,'roi'))
  rmappdata(handles.figure1,'roi');

end

if (isappdata(handles.figure1,'tissue_BW'))
  rmappdata(handles.figure1,'tissue_BW');
rmappdata(handles.figure1,'tissue2_BW');
end

h1= findobj(hAxes,'Tag','roi_pos');
h2= findobj(hAxes,'Tag','inner_pos');
h3= findobj(hAxes,'Tag','outer_pos');

    delete(h1);
    if ~isempty(h2)
     delete(h2);
    end
    if ~isempty(h3)
     delete(h3);
    end
set(handles.clear,'Enable','off');
%set(handles.save,'Enable','off');
set(handles.saveas,'Enable','off');
%set(handles.autoroi,'Enable','off');
%set(handles.autoroi2,'Enable','off');



% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

o_image = getappdata(handles.figure1,'o_image');
cm=getappdata(handles.figure1,'colormap');

rect_pos_new = getappdata(handles.figure1,'rect_pos');
hAx = getappdata(handles.figure1,'axes');

c_dir = pwd;


[fname,dname]=uigetfile('*.mat','Select the roi file','ROIs.mat');
cd(c_dir);

if ischar(fname)
 
     d=whos('-file',fullfile(dname,fname),'roi_pos');
 if ~isempty(d)   
   load(fullfile(dname,fname),'roi','inner_pos','outer_pos','roi_pos','rect_pos');

   
  
   clear_Callback(hObject, eventdata, handles);

  
      roi_pos(:,1) = roi_pos(:,1)+rect_pos(1)-rect_pos_new(1);
      roi_pos(:,2) = roi_pos(:,2)+rect_pos(2)-rect_pos_new(2);
      h=impoly(hAx,roi_pos);
      set(h,'Tag','roi_pos');
      api = iptgetapi(h);
      api.setColor('green');
                
 

if ~isempty(inner_pos)
  inner_pos(:,1) = inner_pos(:,1)+rect_pos(1)-rect_pos_new(1);
  inner_pos(:,2) = inner_pos(:,2)+rect_pos(2)-rect_pos_new(2);
  h=impoly(hAx,inner_pos);
  set(h,'Tag','inner_pos');
  api = iptgetapi(h);
  api.setColor('red');
end
if ~isempty(outer_pos)
  outer_pos(:,1) = outer_pos(:,1)+rect_pos(1)-rect_pos_new(1);
  outer_pos(:,2) = outer_pos(:,2)+rect_pos(2)-rect_pos_new(2);
  h=impoly(hAx,outer_pos);
  set(h,'Tag','outer_pos');
  api = iptgetapi(h);
  api.setColor('red');
end
    
   if    ~isempty(inner_pos) && ~isempty(outer_pos)
      set(handles.autoroi,'Enable','on');
      set(handles.autoroi2,'Enable','on');
   end
      set(handles.saveas,'Enable','on');
        set(handles.clear,'Enable','on');
%      set(handles.save,'Enable','on');

 else
    load(fullfile(dname,fname),'roi');
    roi=image_transform(roi);
    m=bwmorph((roi>0),'remove');
   m=m(rect_pos_new(1):rect_pos_new(1)+size(o_image,1)-1,rect_pos_new(2):rect_pos_new(2)+size(o_image,2)-1);
    [d,cm]=combine_over_under(o_image,m,cm,[1,0,0],m);
    
    imshow(d,cm);
 end

        
end
o_hAx=getappdata(handles.figure1,'o_hAx');
    slice=getappdata(o_hAx,'slice');
    title(sprintf('Slice %d',slice));


% --- Executes on button press in autoroi.
function autoroi_Callback(hObject, eventdata, handles)
% hObject    handle to autoroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.autoroi,'Enable','off');

c_image = getappdata(handles.figure1,'c_image');
 hAx = getappdata(handles.figure1,'axes');

h = findobj(hAx,'Tag','roi_pos');
roi = getappdata(handles.figure1,'roi');

if ~isempty(h)
  api = iptgetapi(h);
  roi_pos = api.getPosition();
  roi=roipoly(c_image, roi_pos(:,1)', roi_pos(:,2)');
  delete(h);
elseif isempty(roi) 
    msgbox('No initial roi defined');
    return;
   
end

h = findobj(hAx,'Tag','inner_pos');
api = iptgetapi(h);
inner_pos = api.getPosition();
inner = roipoly(c_image, inner_pos(:,1)', inner_pos(:,2)');
  
h = findobj(hAx,'Tag','outer_pos');
api = iptgetapi(h);
outer_pos = api.getPosition();
outer = roipoly(c_image, outer_pos(:,1)', outer_pos(:,2)');
                    
h=findobj(handles.figure1,'Tag','darkroi');
dark_roi=get(h,'Value');
%fin_ROI = mdetect3(roi,inner,outer,c_image,dark_roi);    
fin_ROI=mdetect_kmeans(inner,outer,c_image,dark_roi);
add_contour(handles,fin_ROI);
setappdata(handles.figure1,'roi',fin_ROI);
set(handles.autoroi,'Enable','on');

                

% --- Executes on button press in saveas.
function saveas_Callback(hObject, eventdata, handles)
% hObject    handle to saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
                hAx = getappdata(handles.figure1,'axes');
                c_image = getappdata(handles.figure1,'o_image');
                rect_pos = getappdata(handles.figure1,'rect_pos');
                      h = findobj(hAx,'Tag','roi_pos');
                     if h
                        
                        h = findobj(hAx,'Tag','roi_pos');
                        api = iptgetapi(h);
                        roi_pos = api.getPosition();
                        roi = roipoly(c_image, roi_pos(:,1)', roi_pos(:,2)');
                     elseif isappdata(handles.figure1,'roi')
                             
                             roi = getappdata(handles.figure1,'roi');
                             roi_pos = [];
                     else
                         disp('No roi to save');
                     
                     end
                     
                         h = findobj(hAx,'Tag','inner_pos'); 
                         if ~isempty(h)
                           api = iptgetapi(h);
                           inner_pos = api.getPosition();
                           inner = roipoly(c_image, inner_pos(:,1)', inner_pos(:,2)');
                         else
                             inner=[];
                             inner_pos=[];
                         end
                        h = findobj(hAx,'Tag','outer_pos');
                        if ~isempty(h)
                          api = iptgetapi(h);
                          outer_pos = api.getPosition();
                          outer = roipoly(c_image, outer_pos(:,1)', outer_pos(:,2)');
                        else
                           outer = [];
                           outer_pos=[];
                        end
                        roi=image_transform_back(roi);               
                    
hAx=getappdata(handles.figure1,'o_hAx');
slice=getappdata(hAx,'slice');
                   [filename,pathname] = uiputfile('*.mat','Specify the file name','ROIs_temp.mat');
                   if ~isempty(filename)
                       prefix=strtok(filename,'.');
                       prefix=strtok(prefix,'X');
                       filename=[prefix,'X',num2str(slice)];
                      save(fullfile(pathname,filename),'roi','inner','outer','inner_pos','outer_pos','roi_pos','rect_pos');
                      fprintf('ROI saved in %s\n', filename); 
                   
                   end


function  add_contour(varargin)
        
handles = varargin{1};
hAx = getappdata(handles.figure1,'axes');

img = getappdata(handles.figure1,'o_image');
scale=getappdata(handles.figure1,'scale');
img = floor((img-scale(1))*100/diff(scale));
img(img<=2) = 2;

cm=gray(100);
cm(1,:)=[0,1,0];
colormap(hAx,cm);
for i=2:nargin
     
      gap_BW=bwmorph(varargin{i},'remove');
     img(gap_BW>0) = 1;     
 %   img(varargin{i}>0)=1; 
end

h_img = findobj(hAx,'Type','image');
set(h_img,'CData',img);
caxis([1,100]);

function [img,cm]= add_roi_contour(img,cm,handles)
        
roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       h=findobj(handles.figure1,'Tag','roi_pos');
       if isempty(h)
        return;
       end
       api = iptgetapi(h);
       roi_pos = api.getPosition();
       
       
       roi = roipoly(img, roi_pos(:,1)', roi_pos(:,2)');
   end
   
gap_BW=bwmorph(roi,'remove');
img(img>size(cm,1))=size(cm,1);
img(gap_BW>0)=size(cm,1)+1;

cm=cat(1,cm,[0,1,0]);

% --- Executes on button press in autoroi2.
function autoroi2_Callback(hObject, eventdata, handles)
% hObject    handle to autoroi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.autoroi,'Enable','off');

c_image = getappdata(handles.figure1,'c_image');
 hAx = getappdata(handles.figure1,'axes');

h = findobj(hAx,'Tag','roi_pos');
roi = getappdata(handles.figure1,'roi');

if ~isempty(h)
  api = iptgetapi(h);
  roi_pos = api.getPosition();
  roi=roipoly(c_image, roi_pos(:,1)', roi_pos(:,2)');
  delete(h);
elseif isempty(roi) 
    msgbox('No initial roi defined');
    return;
   
end

h = findobj(hAx,'Tag','inner_pos');
api = iptgetapi(h);
inner_pos = api.getPosition();
inner = roipoly(c_image, inner_pos(:,1)', inner_pos(:,2)');
  
h = findobj(hAx,'Tag','outer_pos');
api = iptgetapi(h);
outer_pos = api.getPosition();
outer = roipoly(c_image, outer_pos(:,1)', outer_pos(:,2)');
                    
h=findobj(handles.figure1,'Tag','darkroi');
dark_roi=get(h,'Value');
fin_ROI = mdetect3(roi,inner,outer,c_image,dark_roi);    
%fin_ROI=mdetect_kmeans(inner,outer,c_image,dark_roi);
add_contour(handles,fin_ROI);
setappdata(handles.figure1,'roi',fin_ROI);
set(handles.autoroi,'Enable','on');



% --- Executes on button press in Up.
function Up_Callback(hObject, eventdata, handles)
% hObject    handle to Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       error('no roi found');
   end
  
     roi=circshift(roi,[-1,0]);
     setappdata(handles.figure1,'roi',roi);
      
     add_contour(handles,roi);
     

% --- Executes on button press in Left.
function Left_Callback(hObject, eventdata, handles)
% hObject    handle to Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       error('no roi found');
   end
  
     roi=circshift(roi,[0,-1]);
     setappdata(handles.figure1,'roi',roi);
      
     add_contour(handles,roi);
     

% --- Executes on button press in Down.
function Down_Callback(hObject, eventdata, handles)
% hObject    handle to Down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       error('no roi found');
   end
  
     roi=circshift(roi,[1,0]);
     setappdata(handles.figure1,'roi',roi);
      
     add_contour(handles,roi);

% --- Executes on button press in Right.
function Right_Callback(hObject, eventdata, handles)
% hObject    handle to Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       error('no roi found');
   end
  
     roi=circshift(roi,[0,1]);
     setappdata(handles.figure1,'roi',roi);
      
     add_contour(handles,roi);

% --- Executes on button press in fliplr.
function fliplr_Callback(hObject, eventdata, handles)
% hObject    handle to fliplr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi=getappdata(handles.figure1,'roi');

  
       h=findobj(handles.figure1,'Tag','roi_pos');
       if isempty(h)
        error('no roi found');
       end
       api = iptgetapi(h);
       roi_pos = api.getPosition();
       
       o_image = getappdata(handles.figure1,'o_image');
   
    % roi=fliplr(roi);
    roi_pos(:,1)=size(o_image,2)-roi_pos(:,1);
    api.setPosition(roi_pos);
    
    % setappdata(handles.figure1,'roi',roi);
    
   %  h2=impoly(gca,roi_pos);
   %  add_contour(handles,roi);



% --- Executes on button press in Draw1.
function Draw1_Callback(hObject, eventdata, handles)
% hObject    handle to Draw1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hAx = getappdata(handles.figure1,'axes');


clear_Callback(hObject, eventdata, handles);

title('Draw ROI only');
 
                               
                h = impoly(hAx,[]);
                api = iptgetapi(h);
                api.setColor('green');
                set(h,'Tag','roi_pos');
                
set(handles.clear,'Enable','on'); 
%set(handles.autoroi,'Enable','off');
%set(handles.autoroi2,'Enable','off');
%set(handles.save,'Enable','on');
set(handles.saveas,'Enable','on');


o_hAx=getappdata(handles.figure1,'o_hAx');
    slice=getappdata(o_hAx,'slice');
    title(sprintf('Slice %d',slice));



% --- Executes on button press in rmouterlier.
function rmouterlier_Callback(hObject, eventdata, handles)
% hObject    handle to rmouterlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi=getappdata(handles.figure1,'roi');
c_image = getappdata(handles.figure1,'c_image');

   if isempty(roi)
       h=findobj(handles.figure1,'Tag','roi_pos');
       if isempty(h)
        error('no roi found');
       end
       api = iptgetapi(h);
       roi_pos = api.getPosition();
       
       c_image = getappdata(handles.figure1,'c_image');
       delete(h);
       roi = roipoly(c_image, roi_pos(:,1)', roi_pos(:,2)');
   end
  
   h=findobj(handles.figure1,'Tag','outl_thr');
  
   value=get(h,'String');
   outl_thr =  str2num(value);
   old_roi=roi;
   while 1
       
    mn=mean(c_image(roi>0));
    sd=std(c_image(roi>0));
    roi=roi & c_image>mn-sd*outl_thr & c_image<mn+sd*outl_thr;
    roi = largest_cluster(roi);
    
    
    if length(find(roi>0)) == length(find(old_roi>0))
        break;
    end
    old_roi=roi;
   end
    
    
   setappdata(handles.figure1,'roi',roi);
      
   add_contour(handles,roi);
     

function outl_thr_Callback(hObject, eventdata, handles)
% hObject    handle to outl_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outl_thr as text
%        str2double(get(hObject,'String')) returns contents of outl_thr as a double


% --- Executes during object creation, after setting all properties.
function outl_thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outl_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in darkroi.
function darkroi_Callback(hObject, eventdata, handles)
% hObject    handle to darkroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of darkroi




% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


roi=getappdata(handles.figure1,'roi');

   if isempty(roi)
       h=findobj(handles.figure1,'Tag','roi_pos');
       if isempty(h)
        error('no roi found');
       end
       api = iptgetapi(h);
       roi_pos = api.getPosition();
       
       c_image = getappdata(handles.figure1,'c_image');
       
       roi = roipoly(c_image, roi_pos(:,1)', roi_pos(:,2)');
   end
  
   
hAx=getappdata(handles.figure1,'o_hAx');
  d=get(hAx,'UserData');
  
  y=mean_roi(d,roi);
  figure;
  plot(y);
  

  
   
   



% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Enable','off');
roi=getappdata(handles.figure1,'roi');

hAxes = getappdata(handles.figure1,'axes');
%{
c_image = getappdata(handles.figure1,'c_image');

scale = getappdata(handles.figure1,'scale');
ch = findobj(hAxes,'Type','image');
colormap(gray(200));
set(ch,'CData',c_image);
caxis(scale);
%}
    c_image=getappdata(handles.figure1,'o_image');

    nd=get(handles.edit3,'String');
    nd=str2double(nd);
    cm=getappdata(handles.figure1,'colormap');
    
    [c_image,cm]=add_roi_contour(c_image,cm,handles);
    
    
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
    imwrite(y,get(handles.edit4,'String'),'TIFF');
    fprintf('%s saved\n',get(handles.edit4,'String'));

set(hObject,'Enable','on');

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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


