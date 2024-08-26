function rp=disp_corr(x,y,type)

x=double(x);
y=double(y);
if exist('type','var')
[r,p]=corr(x(:),y(:),'Type',type);
else
    type='Pearson';
[r,p]=corr(x(:),y(:));
  
end
fprintf('%s coef (n=%d) = %3.2e; p = %3.2e\n',type,length(x(:)),r,p);

 rp=[r,p]; 