function swapxy(h)
b=findobj(h,'Tag','swapxy');

cb=get(b,'Callback');

cb(b,'');
showfull(h);
