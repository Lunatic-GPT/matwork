function res=cell2stdarray(c,flag,dim)

if ~exist('flag','var')
    flag='omitnan';
end

res=[];
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
            if isempty(c{i,j,k})
                continue;
            end
            if exist('dim','var')
                res(i,j,k)=std(c{i,j,k},[],dim,flag);
            else
                res(i,j,k)=std(c{i,j,k},flag);
                
            end
        end
    end
    
end

res=squeeze(res);