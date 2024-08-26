function showfull(h)
b=findobj(h,'Tag','showFull');

cb=get(b,'Callback');

cb(b,'');