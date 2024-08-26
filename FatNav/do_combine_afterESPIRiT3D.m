function do_combine_afterESPIRiT3D(fname_raw,nmaps,s)

prefix=strtok(fname_raw,'.');
s=str2num(s);
prefix2=sprintf('%s_ESPIRiT_%dmaps',prefix,nmaps);

wait_for_files(prefix2,s);

d=[];
for i=1:length(s)
    fname=sprintf('%s_%d_*.mat',prefix2,s(i));
    
    fname=dir(fname);
    
    if length(fname)>1
        disp(fname(1).name);
        error('more than one file found');

    end
    
    fname=fname.name;
    fprintf('Loading %s\n',fname);
    dtmp=load(fname);
    d=cat(2,d,dtmp.res);
    
end

fname=sprintf('%s_ESPIRiT_%dmaps.mat',prefix,nmaps);
save(fname,'d','-v7.3');


function wait_for_files(prefix,ind1,nsecond)
found =false;

if ~exist('nsecond','var')
    nsecond = 60;
end

while ~found
   
    found=true;
    for i=1:length(ind1)
        
        name=sprintf('%s_%d_*.mat',prefix,ind1(i));
        dir_str=dir(name);
       
        if isempty(dir_str) || dir_str(1).bytes<1000
            found = false;
            fprintf('%s not found\n',name);
            break;
        end
        
    end
    if ~found
     pause(nsecond);
    end
end