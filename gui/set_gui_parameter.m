function set_gui_parameter(h,tag,val)

d=findobj(h,'Tag',tag);

if ~isempty(d)
    if isa(val,'char')
      set(d,'String',val);
    end
else
    fprintf('Object %s not found\n',tag);
end