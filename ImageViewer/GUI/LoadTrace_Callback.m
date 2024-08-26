% --- Executes on button press in LoadTrace.
function LoadTrace_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


use_pwd=get(handles.pwd,'Value');

tracedirpath=fullfile(tood,'Projects_ongoing','DMV','tracedirpath.mat');
dataDir = pwd;
if ~use_pwd
    if exist(tracedirpath,'file')
      load(tracedirpath);
    end    
end

exts={'*.swc'};

[fname,dataDir]=uigetfile(exts,'Select the roi file',dataDir);

if ~dataDir
    return;
end

save(tracedirpath,'dataDir');

ftrace=fullfile(dataDir,fname);

roi=getappdata(gcf,'uroi');
if isempty(roi)
    data=getappdata(gcf,'udata');
    roi=0*data;
end

dim=[0.4297 0.4297 0.4000];
lpos=read_swc(ftrace,dim);
    lpos{1}(:,1)=417-lpos{1}(:,1);
    lpos{1}(:,2)=513-lpos{1}(:,2);
      
      
p1=round(lpos{1}(1,:));
p2=round(lpos{1}(end,:));

pos=round(connect_pos(lpos{1}));
c=constants;
roi=setv_points(roi,pos,c.vmask);
setappdata(gcf,'TracePath',{p1,pos(2:end,:)});
setappdata(gcf,'TracePoints',[p1;p2]);

setappdata(gcf,'TraceMask',{p1,pos(2:end,:)});
setappdata(gcf,'uroi',roi);

thresholdSlider_Callback(hObject, eventdata, handles);
%showImages(handles);


function pos2=connect_pos(pos)

pos2=[];
for i=1:size(pos,1)-1
    
    n=pos(i+1,:)-pos(i,:);
    
    npix=ceil(sqrt(sum(n.^2)));
    
    for j=0:npix-1
        tmp=round(pos(i,:)+j*n/npix);
        if isempty(pos2) || any(tmp~=pos2(end,:))
            pos2(end+1,:) = tmp;
        end
    end
    
    
end

if  any(pos(end,:)~=pos2(end,:))
    pos2(end+1,:) = pos(end,:);
end