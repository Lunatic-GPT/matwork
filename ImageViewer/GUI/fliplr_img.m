function fliplr_img(h)
b=findobj(h,'Tag','fliplr_image');

cb=get(b,'Callback');

cb(b,'');
showfull(h);
