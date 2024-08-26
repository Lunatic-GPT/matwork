function tstr = fname2title(fname)

[d,f]=fileparts(fname);

[s,t,e] = fnameparts(f);

tstr = sprintf('Slice %d',s);

if ~isempty(t)
    tstr = sprintf('%s, t=%d',tstr,t);
end

if ~isempty(e)
   
    tstr = sprintf('%s (%s)',tstr,e);
    
end

tstr = strrep(tstr,'_','\_');
title(tstr);