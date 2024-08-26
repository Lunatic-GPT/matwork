
function setValue_ui(handle,tag,value)

a=findobj(handle,'Tag',tag);
set(a,'Value',value);