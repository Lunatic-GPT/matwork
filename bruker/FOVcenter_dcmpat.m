function FOVcenter_dcmpat(dname,sl_first)
%dcm2analyze(dname,pattern,nslices,nt)


tic;

if ~exist('sl_first','var')
    sl_first=true;
end

tr=[];
img=[];
pos=[];
dir_str=dir([dname,'/*.dcm']);

bs=[];
es=[];
for i=1:length(dir_str(1).name)
   bs_end=0;
    for j=2:length(dir_str)
       
        if (dir_str(j).name(i)~=dir_str(1).name(i))
           bs_end=1;
            break;
        end 
    end
    
    
        
    if bs_end==1
        break;
    end
    bs(end+1)=dir_str(1).name(i);
    
end
bs=sprintf('%c',bs);

for i=1:length(dir_str(1).name)
   es_end=0;
    for j=2:length(dir_str)
       
        if (dir_str(j).name(end-i+1)~=dir_str(1).name(end-i+1))
           es_end=1;
            break;
        end 
    end
    
    
        
    if es_end==1
        break;
    end
    es(end+1)=dir_str(1).name(end-i+1);
    
end
es=sprintf('%c',es);
es=fliplr(es);

samelen=1;
for j=2:length(dir_str)
        if length(dir_str(j).name)~=length(dir_str(1).name)
            samelen=0;
            break;
        end
end
    
   if samelen
       pat=sprintf('%s%%0%dd%s',bs,length(dir_str(1).name)-length(bs)-length(es),es);
   else
      pat=sprintf('%s%%d%s',bs,es);
   end
   
pos=[];
tr=[];
zerobase=false;
if exist(fullfile(dname,sprintf(pat,0)),'file')
    zerobase=true;
end
for i=1:length(dir_str)
   if mod(i,10)==0
       disp(i);
   end
    if zerobase
        j=i-1;
    else
        j=i;
    end
    
   if ~exist('ns','var')
    h=dicominfo(fullfile(dname,sprintf(pat,j)));
        
    if sl_first
        tmp=h.SliceLocation;
        if ~any(tmp==pos)
         pos(i)=tmp;
        else
            ns=length(pos);
            break;
        end
        
    else
        
    end
    
   end

end

pos2=unique(pos);
ns=length(pos2);


fprintf('offset = %f; nslice = %d; thickness = %f\n',mean(pos2),ns,abs(pos2(2)-pos2(1)));

    disp([mfilename ' finish in ', num2str(toc), ' s']);