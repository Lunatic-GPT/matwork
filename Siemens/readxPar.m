function res=readxPar(fname,par)
%res=readbPar(fname,par,isnum)

fid=fopen(fname,'r');
res=[];
while ~feof(fid)
    
  b=fgetl(fid);
  %{
  if b==-1
      fclose(fid);
      if isempty(res3)
        error([par, ' not found']);
      
      else
       return;    
      end
  end
  %}

  ind=strfind(b,par);
  if ~isempty(ind) 
    
      
    b=fgetl(fid);
     while isempty(strfind(b,'{'))
         b=fgetl(fid);
     end
     
    
     line=1;
     b=fgetl(fid);
     
     while isempty(strfind(b,'}'))
         
         res{line}=strtrim(b);
         line=line+1;
         b=fgetl(fid);
     end
     
     
  else
      continue;
  end

end
if isempty(res)
    disp([par, ' not found']);
end

fclose(fid);

