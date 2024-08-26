function swapxz(h)
b=findobj(h,'Tag','swapxz');

cb=get(b,'Callback');

cb(b,'');
showfull(h);
