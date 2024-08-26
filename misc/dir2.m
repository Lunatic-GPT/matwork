function str=dir2(pattern,dir_only)
%dir_only: 0 - everything; 1 - dir only; 2 - files only
if ~exist('dir_only','var')
    dir_only=false;
end

if ~exist('pattern','var')
    pattern='*';
end

if pattern(end)==':' %disk number
   pattern(end+1)='\';
end

str=dir(pattern);

ind_exc=[];

for i=1:length(str)
  if strcmp(str(i).name,'.') || strcmp(str(i).name,'..') || strcmp(str(i).name,'$RECYCLE.BIN') 
     ind_exc(end+1)=i;
  end
  
end

str(ind_exc)=[];
%{
if strcmp(pattern,'*')
    str=str(3:end); %skip '.' and '..';

elseif length(pattern)>1 && (strcmp(pattern(end-1:end),'\*') || strcmp(pattern(end-1:end),'/*'))
    str=str(3:end);
elseif exist(pattern,'dir')
    str=str(3:end);
end

%}
if dir_only==1
    
    ind_rm=[];
    for i=1:length(str)
       if  ~str(i).isdir
        ind_rm(end+1)=i;
           
       end
    end
    
    str(ind_rm)=[];
elseif dir_only==2
    ind_rm=[];
    for i=1:length(str)
       if  str(i).isdir
        ind_rm(end+1)=i;
           
       end
    end
    
    str(ind_rm)=[];
    
end

