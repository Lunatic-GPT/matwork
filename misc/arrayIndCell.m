function res=arrayIndCell(ind,cl)
%% res=arrayIndCell(ind,cl)
% ind: indices obtained from cell2array
% cl: the cell
% res: length(ind)*2; contains the cell and array indices 


start=0;
for i=2:length(cl)
    
    start(i)=start(i-1)+length(cl{i-1});
end


res=zeros(length(ind),2);  % cell index; element index within the cell

for j=1:length(ind)
    
    tmp=find(ind(j)>start);
    res(j,1)= max(tmp);
    res(j,2)= ind(j)-start(max(tmp));
    
end
