function res=maxlen_cell(c)

res=length(c{1});
for i=2:length(c)
   
    if res<length(c{i})
       
        res=length(c{i});
        
    end
    
    
end