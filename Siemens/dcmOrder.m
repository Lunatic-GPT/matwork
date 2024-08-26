function dcmOrder(path,pattern,timeOrder,doCheck)
% dcmOrder(path,pattern,timeOrder,[doCheck])
% rename the dicom files matched by "pattern" with 0001.ima, 0002.ima ...
% in increasing order of acquisition time.
%if doCheck is true (default), will only report whether the file names are already 
% in the right order but will not rename the files.
% timeOrder: 'z' (slice) or 't' (time)
if ~exist('doCheck','var')
 doCheck = true;
end
cur_d = pwd;
cd(path);
flist = dir2(pattern);
warning off;

for i = 1:length(flist)
in = dicominfo(flist(i).name);
if timeOrder=='z'
 time(i)=in.SliceLocation;
 
elseif timeOrder=='t'
 time(i)=str2double(in.ContentTime);
end    
 
te(i)=in.EchoTime;
end

ute=unique(te);
for ite=1:length(ute)
    fprintf('TE = %f ms\n',te(ite));
    if  timeOrder=='z'
        [stime,ind] = sort(time(te==ute(ite)),'descend');
    elseif timeOrder=='t'
        [stime,ind] = sort(time(te==ute(ite)),'ascend');
    else
        error('unknown type');
    end
    
    disp(['In ' path]);
    
    if isempty(find(diff(ind)<=0,1)) || isempty(find(diff(ind)>=0,1))
        disp('Files names alreay in the right order');
        
    ordered=true;
    else
        disp('Files names not in the right order.');
        disp('File time order:');
        disp(ind);
        
        ordered=false;
    end
    
        if ~doCheck && (~ordered||length(ute)>1)
            disp('Renaming the files ...');
            
            flist_tmp=flist(te==ute(ite));
            newlist = flist_tmp(ind);
            
            ro_ind=strfind(newlist(1).name,'ro');
            
            prefix=['r'*ones(1,ro_ind),'ro'];
            
            for i=1:length(flist_tmp)
                if(mod(i,10)==0)
                    fprintf('%d \n',i);
                end
                
                if length(ute)==1
                  movefile(newlist(i).name,sprintf('%s%04d.dcm',prefix,i));
                else
                  if i==1
                      mkdir(sprintf('TE%dms',round(ute(ite))));
                  end
                  movefile(newlist(i).name,sprintf('TE%dms/%s%04d.dcm',round(ute(ite)),prefix,i));
                  
                      
                end
                    
            end
            
            disp('Done');
        end
        
    end

    


cd(cur_d);

%i1= 44;
%i2 = 14;
%i3= 5;

%convert the above index into two index
%for orientation in rpi
%to3d = zeros(144,1);
%time = zeros(144,1);
%col = mod(i3-1,6)*64+i1;
%row = floor((i3-1)/6)*64+(65-i2);