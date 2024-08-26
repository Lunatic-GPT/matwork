function res=cellCombine(c,dim_cell,dim_array)
% dim_cell: the cell dimension for combining the data; that dimension will
% have only one element after combining
% dim_array: the array dimennsion for combining the data

sz=size(c);
sz(dim_cell)=1;

res=cell(sz);
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
            
            if dim_cell==1
             res{1,j,k}=cat(dim_array,res{1,j,k},c{i,j,k});
            elseif dim_cell==2
             res{i,1,k}=cat(dim_array,res{i,1,k},c{i,j,k});
            else
             res{i,j,1}=cat(dim_array,res{i,j,1},c{i,j,k});
            end
            
            
        end
    end
    
end

res=squeeze(res);