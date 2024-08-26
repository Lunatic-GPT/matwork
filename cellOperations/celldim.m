function dim=celldim(c)

dim=zeros(size(c));
for i=1:size(c,1)
    for j=1:size(c,2)
        dim(i,j)=length(c{i,j});
    end
end