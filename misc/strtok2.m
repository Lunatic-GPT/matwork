function [pre,suf]=strtok2(str,tok)
% backward strtok

str=str(end:-1:1);
[suf,pre]=strtok(str,tok);

if ~isempty(pre)
  pre=pre(end:-1:2);
  suf=[tok,suf(end:-1:1)];
else
    [pre,suf]=strtok(str(end:-1:1),tok);   
end



