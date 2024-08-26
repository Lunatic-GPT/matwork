function [img,cm]=draw_image_roi_onIMG(udata,uscale,odata,oscale,uroi,oroi)

% [img,cm]=draw_image_roi(img,scale,varargin)

setappdata(1,'udata',udata);
set_gui_parameter(1,'edit3',num2str(uscale));

setappdata(1,'odata',odata);
setappdata(1,'odata_rs',[]);

set_gui_parameter(1,'edit12',num2str(oscale));


setappdata(1,'uroi',uroi);


setappdata(1,'oroi',oroi);





