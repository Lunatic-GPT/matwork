
function setString_ui(handle,tag,value)

a=findobj(handle,'Tag',tag);
set(a,'String',value);