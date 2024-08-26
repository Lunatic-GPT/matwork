function res=cellLen(c,dim)


res=[];
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
            if ~exist('dim','var')
                res(i,j,k)=length(c{i,j,k}(:));
            else
                res(i,j,k)=size(c{i,j,k},dim);
            end
        end
    end
    
end