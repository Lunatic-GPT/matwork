function res=cell2nel(c,dim)


res=[];
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
            
            if exist('dim','var')
                res(i,j,k)=size(c{i,j,k},dim);
            else
                res(i,j,k)=length(c{i,j,k});
                
            end
        end
    end
    
end

res=squeeze(res);