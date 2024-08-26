function [roi,ind]=get_ac_pc(prefix,draw_load)

a=ri(sprintf('%s.nii.gz',prefix));
ind=size(a,1)/2;


   % b=sum(sum(a,2),3);
    
  %  [~,ind]=max(b);
  if draw_load==1
    roi=[];
    sag=squeeze(a(ind,:,:));
    sag=flip(flip(sag,1),2);   
    get_line_roi(sag,prefix);

  else
      
    fname=sprintf('ac_pc_%sb.mat',filename(prefix));
    
    while ~exist(fname,'file')
        pause(0.5);
    end
    
    roi=load(fname);
    roi=roi.d;
    
    roi=flip(flip(roi,1),2);

  end
 
 
function get_line_roi(sag,prefix)

fig=img;
handles=guidata(fig);
set(handles.roiShapes,'Value',2);
set(handles.edit3,'String','0 110');
set(handles.showFull,'Value',1);
setappdata(gcf,'udata',sag);
setappdata(gcf,'uroifile',sprintf('ac_pc_%s.mat',filename(prefix)));
showImages(handles);




