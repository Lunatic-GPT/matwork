function m=fftmat(n,inv)
%  m=fftmat(n,inv)

m=zeros(n,n);
for i=1:n
    for j=1:n
        
        if ~exist('inv','var') || ~inv
            m(i,j) = exp(-1i*2*pi*(i-1)*(j-1)/n);
        else
            m(i,j) = exp(1i*2*pi*(i-1)*(j-1)/n);
        end
    end
end
m=m/sqrt(n);