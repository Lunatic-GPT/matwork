function res=histc2(d,bin)

if isempty(d)
    res=zeros(1,length(bin));
else
    
    res=histc(d,bin);
end