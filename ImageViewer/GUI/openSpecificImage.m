 function openSpecificImage(scale)
    
 
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
                
                
                xlm = get(hAx,'XLim');
                ylm = get(hAx,'YLim');
                xlm(1)=ceil(xlm(1));
                xlm(2)=floor(xlm(2));
                ylm(1)=ceil(ylm(1));
                ylm(2)=floor(ylm(2));
                rect_pos(1)= ceil(xlm(1));
                rect_pos(2)= ceil(ylm(1));
                rect_pos(3)=floor(xlm(2))-ceil(xlm(1));
                rect_pos(4)=floor(ylm(2))-ceil(ylm(1));
                 
                
             
            drawROIs(o_image,rect_pos,hAx,'',scale);
            
            
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