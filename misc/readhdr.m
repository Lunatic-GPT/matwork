function res2=readhdr(fname,par)

fid=fopen(fname,'r');
while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      error([par, ' not found']);
  end
  
  ind=strfind(lower(b),lower(par));
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

end

i2=strfind(res,':=');

res2=res(i2+2:end);

while 1
    if ~isempty(res2) && res2(1)==' '
        res2=res2(2:end);
    else
        break;
    end
end

if ~isnan(str2num(res2))
    res2=str2num(res2);
end
fclose(fid);

            