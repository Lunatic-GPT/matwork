function y=repmat_match_cell2array(d,c)

y=[];
for i=1:length(c)

    n=length(c{i});
    y(end+1:end+n)=d(i);
end
    
    