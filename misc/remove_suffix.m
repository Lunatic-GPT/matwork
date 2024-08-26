function name=remove_suffix(name,suf)

n=length(suf);
if length(name)>=length(suf) && strcmp(name(end-n+1:end),suf)
    name=name(1:end-n);
end

    