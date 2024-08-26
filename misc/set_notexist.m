function params=set_notexist(params,field,value)

if ~isfield(params,field)
    params=setfield(params,field,value);
end

    
