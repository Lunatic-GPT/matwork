function res = times(a,b)

if isa(a,'hdrdcp') == 0
    error('In  A.*B only A can be hdrdcp operator');
end

if ndims(b)==3
    b=reshape(b,[size(b,1),size(b,2),1,size(b,3)]);
end


if ndims(b)==2 || ndims(b)==1 
    b=reshape(b,[1,1,1,length(b)]);
end

if a.adjoint
    
    res=zeros(size(b));
        for i=1:size(b,1)
            for j=1:size(b,2)
                for k=1:size(b,3)
                  ts=b(i,j,k,:);
                  ts=squeeze(ts);
                  res(i,j,k,:)=a.xfm*ts;
                end
            end
        end

    
    
else
    
    res=zeros(size(b));
        for i=1:size(b,1)
            for j=1:size(b,2)
                for k=1:size(b,3)
                  ts=b(i,j,k,:);
                  ts=squeeze(ts);
                  res(i,j,k,:)=a.ixfm*ts;
                end
            end
        end
    
    
end


