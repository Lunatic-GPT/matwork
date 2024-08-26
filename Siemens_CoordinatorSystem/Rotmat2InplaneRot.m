function rad=Rotmat2InplaneRot(mat,rad0)
% InplaneDirections2Rotation(mat,deg0)
% rad0: initial solution in radians
% follow dicom convention;
% from [ImageOrientPatient(4:6);ImageOrientPatient(1:3);norm] calculate the
% inplane rotation angle in radian;
% see also NormInplaneRot2Rotmat
%test for T2C2S and S2C2T and deg>0;

norm=mat(:,3)';
norm=norm/sos(norm);

[m,label]=NormInplaneRot2Rotmat(norm,rad0,false);

x=sum(m(:,1).*mat(:,1));
y=sum(m(:,2).*mat(:,1));

rad=atan2(y,x);  
if strcmp(label,'S>C>T')
  rad=-rad;
elseif strcmp(label,'S>T>C')
   rad=-rad; 
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



function res = sos(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end


res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
