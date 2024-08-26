function res=SliceOrient2Norm(rot1,deg1,rot2,deg2)
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


if strcmp(rot,'T2S')
    norm=[-sin(an),0,cos(an)];
elseif strcmp(rot,'S2T')
    norm=[cos(an),0,-sin(an)];
elseif strcmp(rot,'C2T')
    norm=[0,cos(an),-sin(an)];
elseif strcmp(rot,'T2C')
    norm=[0,-sin(an),cos(an)];
elseif strcmp(rot,'S2C')
    norm=[cos(an),-sin(an),0];
elseif strcmp(rot,'C2S')
    norm=[-sin(an),cos(an),0];
else
    error('unknown code');
end

res=m*norm(:);

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
