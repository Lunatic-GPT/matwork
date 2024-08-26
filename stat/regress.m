function [b,F,F2] = regress(x,Y,cov,c,R)
% [b,F,F2] = regress(X,Y,c[,R,cov])
% c is a boolean column vector used to select the coefficients for testing whether including those elements improves fitting.
% R is a matrix added to test R*b = 0;  Can be a cell for multiple tests.
% when R = diag(c), pval=pval2 and F=F2.  R is empty when not used.
% pval is the p value for test contrast specified in c.
% cov is variance-covariance matrix (can be normalized).
% Y is a matrix of nTR*nvox.


if exist('cov','var')
cov=cov*size(cov,1)/trace(cov);
W     = spm_sqrtm(spm_inv(cov));
W     = W.*(abs(W) > 1e-6);
W  = sparse(W);
x=W*x;
Y=W*Y;
end

%x = X;
%b = inv(x'*x)*x'*Y;
%b = (x'*x)\x'*Y;
b=x\Y;

if exist('c','var')&& ~isempty(c)
    
  if length(c)~=size(x,2)
    error('dimension mismatch between c and X');
  end

RSS =sum((Y-x*b).^2,1);

%r2 = 1-RSS/var(Y)/(length(Y)-1);

%% calculate pval and F when c is a row or column vectr

x1 = x(:,~c);
b1 = inv(x1'*x1)*x1'*Y;
RSS1 =sum((Y-x1*b1).^2,1);
n = size(x,1);
p = size(x,2);
%q = length(find(~c))+1;
%F = (RSS1-RSS)/q/(RSS/(n-p));
q=length(find(c));
F = (RSS1-RSS)/q./(RSS/(n-p));
%for i=1:size(Y,2)
% pval(i) = FTest(q,n-p,F(i));
%end

%% calculate pval and F for general con 
end

if exist('R','var')  && ~isempty(R)
  
 if ~iscell(R)   
   q2 = size(R,1);
   Rs=R;
 
   hs = b;
   sig2 = RSS/(n-p);
   F2=zeros(1,size(Y,2));
   
   for j=1:size(Y,2)
     F2(j) = (Rs*hs(:,j))'*inv(Rs*inv(x'*x)*Rs')*Rs*hs(:,j)/q2./sig2(j);
   end
 else
   F2=zeros(length(R),size(Y,2));
   
   for i=1:length(R)
     q2 = size(R{i},1);
     Rs=R{i};
     hs = b;
     sig2 = RSS/(n-p);
     for j=1:size(Y,2)
       F2(i,j) = (Rs*hs(:,j))'*inv(Rs*inv(x'*x)*Rs')*Rs*hs(:,j)/q2./sig2(j);
     end
      % F2(i,:)=diag((Rs*hs(:,j))'*inv(Rs*inv(x'*x)*Rs')*Rs*hs(:,j))/q2./sig2(j);
      % same, but more memory
   end
     
 end
 
end


