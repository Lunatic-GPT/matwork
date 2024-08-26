function [label,deg2,deg1]=Norm2SliceOrient(norm)
% follow dicom convention;
% abs(deg2)>abs(deg1)

norm=norm(:)';
norm=norm/sos(norm);

[deg2,deg1]=C2S2T(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    res=sprintf('C>S%3.2f>T%3.2f\n',deg2,deg1);
    label='C>S>T';
    return;
end

[deg2,deg1]=C2T2S(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    res=sprintf('C>T%3.2f>S%3.2f\n',deg2,deg1);
    label='C>T>S';
    return;
end

[deg2,deg1]=T2S2C(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    res=sprintf('T>S%3.2f>C%3.2f\n',deg2,deg1);
    label='T>S>C';
    return;
end

[deg2,deg1]=T2C2S(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    res=sprintf('T>C%3.2f>S%3.2f\n',deg2,deg1);
    label='T>C>S';
    return;
end

[deg2,deg1]=S2C2T(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
   res=sprintf('S>C%3.2f>T%3.2f\n',deg2,deg1);
   label='S>C>T';
   return;
end

[deg2,deg1]=S2T2C(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    res=sprintf('S>T%3.2f>C%3.2f\n',deg2,deg1);
    label='S>T>C';
    return;
end
%disp(res);

function [a,b]=C2S2T(norm)
b=-asin(norm(3));
a=atan2(-norm(1),norm(2));
a(~any(abs(norm(1:2))>1e-6))=0;
[a,b]=fix_ab(a,b);

function [a,b]=C2T2S(norm)
b=-asin(norm(1));

a=atan2(-norm(3),norm(2));
a(~any(abs(norm(2:3))>1e-6))=0;
[a,b]=fix_ab(a,b);

function [a,b]=T2S2C(norm)
b=-asin(norm(2));
a=atan2(-norm(1),norm(3));
a(~any(abs(norm([1,3]))>1e-6))=0;
[a,b]=fix_ab(a,b);


function [a,b]=T2C2S(norm)
b=-asin(norm(1));
a=atan2(-norm(2),norm(3));
a(~any(abs(norm([2,3]))>1e-6))=0;
[a,b]=fix_ab(a,b);



function [a,b]=S2T2C(norm)
b=-asin(norm(2));
a=atan2(-norm(3),norm(1));
a(~any(abs(norm([1,3]))>1e-6))=0;
[a,b]=fix_ab(a,b);



function [a,b]=S2C2T(norm)
b=-asin(norm(3));
a=atan2(-norm(2),norm(1));
a(~any(abs(norm([1,2]))>1e-6))=0;
[a,b]=fix_ab(a,b);



function [a,b]=fix_ab(a,b)

a=a*180/pi;
b=b*180/pi;

if a>90
    a=a-180;
    b=-b;
elseif a<-90
    a=a+180;
    b=-b;
end

