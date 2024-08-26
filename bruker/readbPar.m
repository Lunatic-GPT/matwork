function res=readbPar(fname,par,isnum)
%res=readbPar(fname,par,isnum)

if ~exist('isnum','var')
    isnum=true;
end

fid=fopen(fname,'r');
par=['##$',par,'='];
while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      error([par(4:end-1), ' not found']);
  end
  
  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

end

if res(1)=='('
    sz=str2num(res(2:end-1));
     if numel(sz)==1
            sz=[1,sz];
     end
        
    if isnum 
     res=[];
     while 1
        
      b=fgetl(fid);
      if isnum
        tmp=str2num(b);
      else
       tmp=strread(b,'%s');  
      end
      res=[res,tmp];
      
      if length(res)==prod(sz) || ~isnum
          break;
      end
      
     end
    res=reshape(res,sz(end:-1:1));
    
    else
         b=fgetl(fid);
       res=strread(b,'%s');  
     
    end
else
    if isnum    
      res=str2double(res);
    end
    
    
end
fclose(fid);

