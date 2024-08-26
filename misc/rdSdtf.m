function d=rdSdtf(fname)
% specify subbriks with the following formats data+orig[x..y] or data+orig[x]. 0 based
[fname,r] = strtok(fname,'[');
if isempty(r)
    d= rdSdt(fname);
else
    r =strtok(r(2:end),']');
    [r,r2] = strtok(r,'.');
    d=rdSdt(fname);
    
    if isempty(r2)
      r=str2double(r);
      d=d(:,:,r+1);
    else
      r=str2double(r);
      r2=str2double(r2(3:end));
      d= d(:,:,:,r+1:r2+1);
    end
    
end
