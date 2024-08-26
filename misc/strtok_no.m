function [res,suf]=strtok_no(str,tok,rep)


res=[tok,str];

for i=1:rep-1
[tmp, res]=strtok(res(2:end),tok);
end

[res,suf]=strtok(res(2:end),tok);
