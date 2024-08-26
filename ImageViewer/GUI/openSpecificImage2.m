 function openSpecificImage2
    
 
h=findall(gcf,'tag','FigureToolBar');
h2=findall(h,'tag','Exploration.DataCursor');
st=get(h2,'State');

 if strcmp(st,'on')
     return;
 end
 type = get(gcf,'SelectionType');
    
    
    switch type
        case 'normal' % left mouse click             
             % d=read_afni_images(filename);
             
                hAx = get(gcbo,'Parent');
                o_image=get(gcbo,'CData');
             % o_image=d(:,:,slice);
              
              %o_image = double(imread(filename));
                
                cm=getappdata(hAx,'colormap');
                
                 
               lr=findobj(gcf,'Tag','edit14');
               ud=findobj(gcf,'Tag','edit15');
              croplr=get(lr,'String');
              cropud=get(ud,'String');
              croplr=str2num(croplr);
                cropud=str2num(cropud);
       
               
               rect_pos(1)=croplr(1)+1;
               rect_pos(2)=cropud(1)+1;
               
            drawROIs(o_image,cm,rect_pos,hAx,'');
            
            
        case 'open'   
            %double click
        case 'extend'
            % shift & left mouse button action
       
             
        case 'alt'
            % alt & left mouse button action    
            
                hAx = get(gcbo,'Parent');
           % im = imread(filename);
             d=get(hAx,'UserData');
             
              figure;plot(squeeze(d));

    end