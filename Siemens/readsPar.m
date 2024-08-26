function res3=readsPar(fname,par)
%res=readbPar(fname,par,isnum)

fid=fopen(fname,'r');
res3=[];
while ~feof(fid)
    
  b=fgetl(fid);
 

  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    res=strrm(res,' ');
    res=strrm(res,'=');

  else
      continue;
  end


res2=str2double(res);

if ~isnan(res2)    
      res3(end+1)=res2;
else
    res3{end+1}=res;
end
end
if isempty(res3) && nargout==0
    disp([par, ' not found']);
end

fclose(fid);

if isempty(res3)
    res3=0;
end

function oldstr = strrm(oldstr,substr)
%
% newstr = strrm(oldstr,substr)
%remove substring substr from the string oldstr.

while 1
 i=strfind(oldstr,substr);
 if isempty(i)
     break;
 end
 oldstr(i:i+length(substr)-1)=[];
end





