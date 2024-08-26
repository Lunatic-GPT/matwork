function res=get_img_handle

h=findobj('Type','figure');
res=[];
for i=1:length(h)
    if strcmp(get(h(i),'FileName'),'C:\Users\zongx\Dropbox\matwork\blat\image_gallery_afni\imageGallery.fig')
      res(end+1)=h(i);
    end
    
end