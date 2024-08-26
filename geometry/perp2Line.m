function res=perp2Line(v)

% v is a 2 dim vector

if v(1)==0
 res=[1,0];    
else
res=[-v(2)/v(1),1];

res=res/sos(res,2);

end