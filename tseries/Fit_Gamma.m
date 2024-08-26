function beta = Fit_Gamma(x,y,ind_pw)
%beta = Fit_DoubleGamma(x,y,ind_pw)
% ind_pw:  index to fit the positive and negative responses separately.

xx = [];
for i=1:length(x)-1
    xx = [xx,linspace(x(i),x(i+1),11)];
end
xx = unique(xx);
    
b= [5,1.1,-min(y),12,0.9];

ymx = max(y);
ymn = min(y);
options = statset('Robust','on','MaxIter',1000);
if exist('ind_pw','var') && ~isempty(ind_pw)
  beta = zeros(2,2);
  b = [5,1.1];
  beta(:,1) = nlinfit(x(1:ind_pw),y(1:ind_pw)/ymx,@Gamma,b);
  
  b = [12,0.9];
  beta(:,2) = nlinfit(x(ind_pw:end),y(ind_pw:end)/ymn,@Gamma,b);
  
  figure;plot(x,y,'o');
  hold on;plot(xx(1:(ind_pw-1)*10),Gamma(beta(:,1),xx(1:(ind_pw-1)*10))*ymx,'-b');
  
  hold on;plot(xx((ind_pw-1)*10:end),Gamma(beta(:,2),xx((ind_pw-1)*10:end))*ymn,'-r');
  
else
    
  beta = nlinfit(x,y/ymx,@Gamma,b,options);
  figure;plot(x,y,'o',xx,Gamma(beta,xx)*ymx,'-');
  
  
end


    
    

function y = Gamma(b,x)

z = b(1)*b(2);
amp1 = z^b(1)*exp(-z/b(2));

y = x(:).^b(1).*exp(-x(:)/b(2))/amp1;



