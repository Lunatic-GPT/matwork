function cov3=cov_matrix(s)
% s is in the format of nch*npt and can be complex

nch=size(s,1);
cov3=zeros(nch,nch);
    for i=1:nch
        for j=1:nch
            if (i>j)
                continue;
            end
            
            x=s(i,:);
            y=s(j,:);
            cov3(i,j)=mean(conj(x-mean(x)).*(y-mean(y)));
            cov3(j,i)=conj(cov3(i,j));
        end
    end
