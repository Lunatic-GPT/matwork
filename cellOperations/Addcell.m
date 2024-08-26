function res=Addcell(c,dim)
% all cell contents should have the same matrix size

if numel(size(c{1}))>9
    error('Matrix dimensions too large');
end

c=cellCombine(c,dim,10);

sz=size(c);

res=cell(sz);
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
             res{i,j,k}=sum(c{i,j,k},10);
        end
    end
end