function tab=loadtable(fname,tname)

if ~exist('tname','var')
    tname='t1';
end
if exist(fullfile(pwd,fname),'file')
    fid=fopen(fname);
elseif   strcmp(computer,'PCWIN64')
   fid=fopen(['z:/homes/xiaopeng/vnmrsys/tablib/',fname],'r');
elseif exist('/data/users/xiaopeng','dir')
    fid=fopen(['/data/users/xiaopeng/newvnmrsys/tablib/',fname]);
else
   fid=fopen(['/home/xiaopeng/vnmrsys/tablib/',fname]);  
end

while 1
    a=fgetl(fid);
    if feof(fid)
        error(['Table ',tname,' not found/']);
    end
    if ~isempty(findstr(a,tname))
        tab=[];
        while 1
            b = fgetl(fid);
            if isempty(str2double(b)) || isnan(str2double(b))
                break;
            end
            tab(end+1)=str2double(b);
            if  feof(fid)
                break;
            end
        end
        break;
    end
    
end

