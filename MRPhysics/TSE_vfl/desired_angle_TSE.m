function [a,err]=desired_angle_TSE(z1,f1,fm1,s)

%z1=-z1;
tmp=z1*z1-(f1-s)*(fm1-s);

if tmp<0
    a=0;
    err=1;
    return;
    %error(sprintf('negative %f ',tmp));
end
y=sqrt(tmp);


a(1)=2*atan2((z1-y),(f1-s));

% s=verify(a,z1,f1,fm1);
% s=verify2(y2,z1,f1,fm1);


if abs(f1-s)<1e-8   
    a(2)=2*atan2((fm1-f1),2*z1);  %should be fm1-f1
    
elseif abs(fm1-s)<1e-8   
    
    a(2)=2*atan2(-2*z1,(fm1-f1));
else


%y2=(-z1+y)/(f1-s);

%a(2)=2*atan(y2);

a(2)=2*atan2((z1+y),(f1-s));
end

if a(2)<0
    a(2)=a(2)+2*pi;
end
%disp([z1,f1,fm1,s,tmp,a(2)]);
err=0;
%s=verify(a,z1,f1,fm1);
%disp(s);


%s=verify2(y2,z1,f1,fm1);
%disp(s);



%%
function s=verify(a,z1,f1,fm1)

s=(1+cos(a))/2*fm1+(1-cos(a))/2*f1+sin(a)*z1;




function s = verify2(tana2,z1,f1,fm1)

cosa2sq=1/(1+tana2^2);

cosa=2*cosa2sq-1;


sina=2*cosa2sq*tana2;

s=(1+cosa)/2*fm1+(1-cosa)/2*f1+sina*z1;