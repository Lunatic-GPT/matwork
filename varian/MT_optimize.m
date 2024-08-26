%n = [40,40,10,10,10];

%MTR=[0.7,0.53,0.37,0.21,0];

n=[10,10];
MTR = [0.9,0];
x = [1-MTR(:),ones(length(MTR),1)];

N=10000;

b = zeros(1,N);
for i=1:N
    ns = randn_white(length(MTR))./sqrt(n(:));
    tmp=inv(x'*x)*x'*ns;   
    b(i)=tmp(2);
end

fprintf('STD per trial is %f\n',std(b)*sqrt(sum(n)));





