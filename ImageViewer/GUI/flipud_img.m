function flipud_img(h)
b=findobj(h,'Tag','flipud_image');

cb=get(b,'Callback');

cb(b,'');
showfull(h);
