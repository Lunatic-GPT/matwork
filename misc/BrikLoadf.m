function [d,info]=BrikLoadf(fname)
% specify subbriks with the following formats data+orig[x..y] or data+orig[x]. 0 based
[fname,r] = strtok(fname,'[');
if isempty(r)
    [d,info]= BrikLoad(fname);
else
    r =strtok(r(2:end),']');
    [r,r2] = strtok(r,'.');
    if isempty(r2) 
      opt.Frames = str2double(r)+1;
    else
      opt.Frames = str2double(r)+1:str2double(r2(3:end))+1;
    end
    [d,info]=BrikLoad(fname);
    d=d(:,:,:,opt.Frames);
end
