function m2=tabc(m,tabfile,tname)
%m2=tabc(m,tabfile,tname)
% tname: table name. (default t1)
% phase encoding should be the second dimension.
if ~exist('tname','var')
    tname = 't1';
end

fid=fopen(tabfile);
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
        
    

m2 = zeros(size(m));
if length(tab)~=size(m2,2)
    error('table dimension mismatch');
end

[i2,ind] = sort(tab);
for i=1:size(m2,2)
   m2(:,i,:,:) = m(:,ind(i),:,:);
end

    
    

