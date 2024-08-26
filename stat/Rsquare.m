function res=Rsquare(res,d)
% tested in R
% res: difference between fit and data
% d: the data
TSS = sum((d-mean(d)).^2);  %total sum of squares
RSS = sum(res.^2);%residual sum of squares

res=1-RSS/TSS;
