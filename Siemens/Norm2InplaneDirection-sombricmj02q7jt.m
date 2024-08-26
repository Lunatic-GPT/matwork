function res=Norm2InplaneDirection(norm)
% follow dicom convention;

norm=norm(:)';
norm=norm/sos(norm);

[deg2,deg1]=C2S2T(norm);

if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    fprintf('C>S%3.2f>T%3.2f\n',deg2,deg1);
    res=SliceOrient2Mat('C2S',deg2,'C2T',deg1);
end

[deg2,deg1]=C2T2S(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)   
    fprintf('C>T%3.2f>S%3.2f\n',deg2,deg1);
end

[deg2,deg1]=T2S2C(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    fprintf('T>S%3.2f>C%3.2f\n',deg2,deg1);
    res=SliceOrient2Mat('T2S',deg2,'T2C',deg1);
end

[deg2,deg1]=T2C2S(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    fprintf('T>C%3.2f>S%3.2f\n',deg2,deg1);
    res=SliceOrient2Mat('T2C',deg2,'T2S',deg1);
    
end

[deg2,deg1]=S2C2T(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)  

    fprintf('S>C%3.2f>T%3.2f\n',deg2,deg1);
end

[deg2,deg1]=S2T2C(norm);
if abs(deg2)<=45 && abs(deg2)>=abs(deg1)    
    fprintf('S>T%3.2f>C%3.2f\n',deg2,deg1);
end

function [deg1,deg2]= round_angle(deg1,deg2)
    %deg2=round(deg2*10)/10;
    
%deg1=round(deg1*10)/10;


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

function res=SliceOrient2Mat(rot1,deg1,rot2,deg2)
% follow dicom convention;
% the output is the slice direction vector perpendicular to the slice;
% deg1 come first in the string and > deg2
% for S2T2C1.1 rot1='S2T' and rot2='S2C' 

if exist('rot2','var')
   rot=rot2;
   an=deg2*pi/180;
   m=rotmat(rot1,deg1);
else
    rot=rot1; 
    an=deg1*pi/180;
    m=eye(3);
end

m2=rotmat(rot,an*180/pi);

%when inplaneRot = 0
if strcmp(rot,'T2S')
    norm=[-sin(an),0,cos(an)];
    order=[1,2,3];
    m2(:,2)=-m2(:,2);
elseif strcmp(rot,'S2T')  
    norm=[cos(an),0,-sin(an)];
    order=[3,2,1];
elseif strcmp(rot,'C2T')
    norm=[0,cos(an),-sin(an)];
    order=[ 1,3,2];
elseif strcmp(rot,'T2C')
    norm=[0,-sin(an),cos(an)];
    order=[2,1,3];
elseif strcmp(rot,'S2C')
    norm=[cos(an),-sin(an),0];
    order=[3,2,1];
elseif strcmp(rot,'C2S')
    norm=[-sin(an),cos(an),0];
    order=[1,3,2];
else
    error('unknown code');
end

%res=m*norm(:);
res=m*m2(:,order);

%res=norm*m;




function m=rotmat(rot,deg)

an=deg*pi/180;
if strcmp(rot,'T2S')
    x=[cos(an),0,sin(an)];
    y=[0,1,0];
    z=[-sin(an),0,cos(an)];
    
elseif strcmp(rot,'S2T')
    
    x=[cos(an),0,-sin(an)];
    y=[0,1,0];
    z=[sin(an),0,cos(an)];
    
elseif strcmp(rot,'C2T')
    x=[1,0,0];
    y=[0,cos(an),-sin(an)];
    z=[0,sin(an),cos(an)];
elseif strcmp(rot,'T2C')
   x=[1,0,0];
   y=[0,cos(an),sin(an)];
   z=[0,-sin(an),cos(an)];
    
elseif strcmp(rot,'S2C')
    
    x=[cos(an),-sin(an),0];
    y=[sin(an),cos(an),0];
    z=[0,0,1];
elseif strcmp(rot,'C2S')
    x=[cos(an),sin(an),0];
    y=[-sin(an),cos(an),0];
    z=[0,0,1];
    
else
    error('unknown code');
end
m=[x',y',z'];

