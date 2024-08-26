function B1map_varFA(dlist,falist)


d=[];
for i=1:length(dlist)
    d(:,:,:,i)=ri(num2str(dlist(i)));
    
end

    
dmax=max(d(:));
    b1=zeros(size(d(:,:,:,1)));
    
for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)
            y=d(i,j,k,:);
            y=y(:);
            
            if max(y)<dmax*0.03
                continue;
            end
            
            x=falist;
            f0=pi/180;
            
            if i==31 && j==75
                disp('');
            end
            
            [b,r] = nlinfit(x(:),y(:),@cosfita,[max(y),f0]);
            
            b1(i,j,k,1) = 90*b(2)/f0; %convert to hertz at the maximum power used.
            
        end
    end
end

save(sprintf('B1map_varFA_%d',dlist(1)),'b1');

