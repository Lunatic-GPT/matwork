function push_button(h,String)
b=findobj(h,'String',String);

cb=get(b,'Callback');

cb(b,'');
showfull(h);
