function res=cell2arrayMax(c,n)

if ~exist('n','var')
    n=1;
end

res=zeros(size(c));
for i=1:size(c,1)
    for j=1:size(c,2)
        for k=1:size(c,3)
            if isempty(c{i,j,k})
                continue;
            end
            %if exist('dim','var')
              %  res(i,j,k,:)=max(c{i,j,k},[],dim);
            %else
   
                tmp=sort(c{i,j,k});
                tmp(isnan(tmp))=[];
                
                if isempty(tmp)
                    continue;
                elseif length(tmp)<n
                    res(i,j,k)=mean(tmp);
                   % continue;
                else
                   res(i,j,k)=mean(tmp(end-n+1:end));
                end
            %end
            
        end
    end
    
end

res=squeeze(res);