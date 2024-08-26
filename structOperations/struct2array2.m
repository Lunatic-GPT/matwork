function res=struct2array2(str,fd)

for i=1:size(str,1)
    for j=1:size(str,2)
      res(i,j)=getfield(str(i,j),fd);
         
        
    end
end