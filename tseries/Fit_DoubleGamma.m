function [beta,t_min,t_max,y_min,y_max] = Fit_DoubleGamma(x,y,b,ind_pw)
%[beta,t_min,t_max,y_min,y_max] = Fit_DoubleGamma(x,y,b,ind_pw)
% b initial values of the fitting parameters.
% ind_pw:  index to fit the positive and negative responses separately.
mx = max(y);
y = y/max(y);

xx = [];
for i=1:length(x)-1
    xx = [xx,linspace(x(i),x(i+1),11)];
end
xx = unique(xx);
yy = zeros(size(xx));   
if ~exist('b','var') || isempty(b) 
  b= [9,0.7,-min(y),7,2,1];
end

 figure;plot(x,y,'o');
options = statset('Robust','on','MaxIter',5000);
if exist('ind_pw','var') && ~isempty(ind_pw)
  beta = zeros(6,2);
  beta(:,1) = nlinfit(x(1:ind_pw),y(1:ind_pw),@doubleGamma,b);
  beta(:,2) = nlinfit(x(ind_pw:end),y(ind_pw:end),@doubleGamma,b);
  figure;plot(x,y,'o');
  yy(1:(ind_pw-1)*10) = doubleGamma(beta(:,1),xx(1:(ind_pw-1)*10));
  yy ((ind_pw-1)*10+1:end) = doubleGamma(beta(:,2),xx((ind_pw-1)*10+1:end));
  hold on;plot(xx,yy,'-b');
  
else
    
  beta = nlinfit(x,y,@doubleGamma,b,options);
  yy = doubleGamma(beta,xx);
 
  hold on;plot(xx,yy,'-');
end

[y_max,t_max] = max(yy);
y_max = y_max*mx;
[y_min,t_min] = min(yy);
y_min = y_min*mx;
t_max = xx(t_max);
t_min = xx(t_min);

function y = doubleGamma(b,x)

z1 = b(1)*b(2);
z2 = b(4)*b(5);

amp1 = z1^b(1)*exp(-z1/b(2));
amp2 = z2^b(4)*exp(-z2/b(5));

%y = b(1)*x(:).^b(2).*exp(-x(:)/b(3)) - b(4)*x(:).^b(5).*exp(-x(:)/b(6));
y = x(:).^b(1).*exp(-x(:)/b(2))/amp1-b(3)*x(:).^b(4).*exp(-x(:)/b(5))/amp2;
y = y*b(6);





