function s=randn_correlated(n,alpha)
% s=randn_correlated(n[,alpha])
% default: alpha = 0.4; correlation decays as alpha^n;

if ~exist('alpha','var')
    alpha=0.4;
end

if alpha>=1  || alpha<0
    error('alpha should be >0 and <1');
end
    
if length(n)==1
    n=[1,n];
end
z=zeros(n);

z(1)=randn_white(1);
for i=2:length(z(:))
    z(i)=alpha*z(i-1)+randn_white(1);
end
  
s=z*sqrt(1-alpha^2);




