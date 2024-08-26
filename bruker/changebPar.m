function res=changebPar(fname,par,isnum,newvalue)
%res=readbPar(fname,par,isnum)

[pathstr,filename]=fileparts(fname);
fid2=fopen(fullfile(pathstr,'temp'),'w');
%disp('create file temp');
if ~exist('isnum','var')
    isnum=true;
end

fid=fopen(fname,'r');
par=['##$',par,'='];
while 1
  b=fgetl(fid);
  
  if b==-1
      fclose(fid);
      fclose(fid2);
      error([par(4:end-1), ' not found']);
  end
  
  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

  fprintf(fid2,'%s\n',b);
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
    if size(newvalue,1)==1 || size(newvalue,2)==1 
     fprintf(fid2,'%s( %d )\n',par,length(newvalue));
    else
     fprintf(fid2,'%s( %d, %d )\n',par,size(newvalue,2),size(newvalue,1));
    end 
    fprintf('setting %s\n',par);
    for i=1:length(newvalue(:))
        fprintf(fid2,'%d ',newvalue(i));
        
        fprintf('%d ',newvalue(i));
    end
      fprintf(fid2,'\n');
      
      fprintf('\n');
    else
         b=fgetl(fid);
       res=strread(b,'%s');  
      
    fprintf(fid2,'%s( %d )\n',par,length(newvalue));
    
    fprintf('setting %s \n',par);
    for i=1:length(newvalue(:))
      fprintf(fid2,'%s ',newvalue{i});
      
      fprintf('%s ',newvalue{i});
    end
      fprintf(fid2,'\n');
      fprintf('\n');
    end
else
  if isnum
    fprintf(fid2,'%s%d\n',par,newvalue);
      fprintf('setting %s%d\n',par,newvalue);
  else
    fprintf(fid2,'%s%s\n',par,newvalue);
    
    fprintf('setting %s%s\n',par,newvalue);
  end
  
end

while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      fclose(fid2);
      if ~exist(fullfile(pathstr,[filename,'.orig']),'file')
        movefile(fname,fullfile(pathstr,[filename,'.orig']),'f');
      end
    %  fprintf('move %s to %s\n',fname,fullfile(pathstr,[filename,'.orig']));
      movefile(fullfile(pathstr,'temp'),fname,'f');
    %  fprintf('move %s to %s\n',fullfile(pathstr,'temp'),fname);
    break;
  else
        fprintf(fid2,'%s\n',b);
  end
  

end




