function [res,label]=NormInplaneRot2Rotmat(norm,rad,apply_inplaneRot)
% [res,label]=NormInplaneRot2Rotmat(norm,rad,apply_inplaneRot)
% follow dicom convention;
% 
% if apply_inplaneRot is false, in plnae rotation will not be applied.
% norm = res(:,3);
% tested for all orientations, but results may not be correct when inplane
% angle is close to 45, 135, 225, or 315.
% under Projects_ongoing\Diabetes\ProtocolOptimization\TestInplaneOrient

%Input: 
% norm: the slice normal direction
% rad: the inplane rotation in radians
% apply_inplaneRot: whether apply_inplaneRot when calculate res;
% if rad is empty or not exist, norm
% and rad will be read from a file with name norm;

% Output:
% res: the rotmat calculated from norm and rad. 
% label: the label for slice orientation

% Note: res should be equal to [ImageOrientPatient(4:6);ImageOrientPatient(1:3);norm]; 
% (same as rotmat from ri) where
% ImageOrientationPatient is same as in dicom header

if ~exist('apply_inplaneRot','var')
    apply_inplaneRot=true;    
end

if ~exist('rad','var') || isempty(rad)
    if length(norm)<=4 || ~strcmp(norm(end-3:end),'.pro')
    a=readdPar(norm,'ImageOrientationPatient');
    fprintf('ImageOrientationPatient in dicom header %s\n',num2str(a(:)'));
    extp(norm);
    fname=[norm,'.pro'];
    else
        fname=norm;
    end
    
    x=readsPar(fname,'sSliceArray.asSlice[0].sNormal.dSag');
    y=readsPar(fname,'sSliceArray.asSlice[0].sNormal.dCor');
    z=readsPar(fname,'sSliceArray.asSlice[0].sNormal.dTra');
    rad=readsPar(fname,'sSliceArray.asSlice[0].dInPlaneRot');
    norm=[x,y,z];
    fprintf('Read Normal and Inplane rotation angle from pro: norm = %s; angle = %f deg\n',num2str(norm(:)'),rad*180/pi);
    
end

deg=rad*180/pi;
norm=norm(:)';
norm=norm/sos(norm);

[label,deg2,deg1]=Norm2SliceOrient(norm);

if strcmp(label,'C>S>T') %tested
    res=SliceOrient2Mat('C2S',deg2,'C2T',deg1);
  [deg,res]=do_swap_neg_m90(deg,res);
    
elseif strcmp(label,'C>T>S')  %tested;
    res=SliceOrient2Mat('C2T',deg2,'C2S',deg1);
    
  [deg,res]=do_swap_neg_m90(deg,res);
    
elseif strcmp(label, 'T>S>C') %tested
     res=SliceOrient2Mat('T2S',deg2,'T2C',deg1);
   [deg,res]=do_swap_neg_0(deg,res);
   
elseif strcmp(label,'T>C>S') %tested
    res=SliceOrient2Mat('T2C',deg2,'T2S',deg1);
   [deg,res]=do_swap_neg_m90(deg,res);
elseif strcmp(label,'S>C>T') %tested
    
    res=SliceOrient2Mat('S2C',deg2,'S2T',deg1);
    if closeto(deg,0) %tested 30 
        res(:,1)=-res(:,1);    %new
    elseif closeto(deg,90) %tested 120
        res(:,[1,2])=res(:,[2,1]); %new
    elseif closeto(deg,-90) %tested -80
        res(:,[1,2])=-res(:,[2,1]);
    else
        res(:,2)=-res(:,2);
    end
    deg=-deg;
else  %not tested S>T>C
    res=SliceOrient2Mat('S2T',deg2,'S2C',deg1);
    if closeto(deg,0) %tested 30 
        res(:,1)=-res(:,1);    %new
    elseif closeto(deg,90) %tested 120
        res(:,[1,2])=res(:,[2,1]); %new
    elseif closeto(deg,-90) %tested -80
        res(:,[1,2])=-res(:,[2,1]);
    else
        res(:,2)=-res(:,2);
    end
    deg=-deg; 
%     if closeto(deg,0)  
%         res(:,1)=-res(:,1);    %new
%     elseif closeto(deg,90) 
%         res(:,[1,2])=res(:,[2,1]); %new
%     elseif closeto(deg,-90)
%         res(:,[1,2])=-res(:,[2,1]);
%     else  %180
%         res(:,2)=-res(:,2);
%     end
    
   % deg=-deg;
    
    
end

if apply_inplaneRot
    [theta,phi]=unitVec2thetaPhi(res(:,3));
    m=transform_matrix_rotation_arb_axis(theta,phi,deg);
    res=m*res;
end

function [deg,res]=do_swap_neg_0(deg,res)
 % no change when angle close to 0 deg
    if closeto(deg,90)  %test 120
        res(:,2)=-res(:,2);
        res(:,[1,2])=res(:,[2,1]); %new 
    elseif closeto(deg,-90) %test -120
        res(:,1)=-res(:,1);
        res(:,[1,2])=res(:,[2,1]); %new
    elseif closeto(deg,180) 
        res(:,[1,2])=-res(:,[1,2]);
    else %test 30 and -30
    end
    deg=-deg;
    
    function [deg,res]=do_swap_neg_m90(deg,res)
        %no change when angle close to -90;
     if closeto(deg,90)  %tested for 120 
        res(:,[1,2])=-res(:,[1,2]);
    elseif closeto(deg,-90) % -46, -60, -45.99, not sure how the threshold is decided.
        
    elseif closeto(deg,0) %tested for 30 and -30, -44, -45.15,-45.16
        res(:,2)=-res(:,2);
        res(:,[1,2])=res(:,[2,1]); %new
    else %180
        res(:,1)=-res(:,1);
        res(:,[1,2])=res(:,[2,1]); %new
     end
    deg=-deg;
    
function res = closeto(deg,ref,range)
% deg is between 0 - 360
% ref also between 0 -360
if ~exist('range','var')
    range=45;
end

d=mod(deg-ref,360);
res= d<=range || (360-d)<=range;
    
    
    
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

%when inplane rotation is 0
if strcmp(rot,'T2S')
    order=[1,2,3];
    m2(:,2)=-m2(:,2); 
    res=m*m2(:,order);
elseif strcmp(rot,'S2T')
    order=[3,2,1];
       res=m*m2(:,order);
elseif strcmp(rot,'C2T')
    order=[ 1,3,2];
       res=m*m2(:,order);
elseif strcmp(rot,'T2C')
    order=[2,1,3];
    res=m*m2;
    
    deg3=-atan2(res(1,2),res(1,1))*180/pi;
    m3=rotmat('S2C',deg3);%somehow need 1.6 deg when deg1 = 15 and deg2 = 6;
    
    res=res*m3(:,order);
elseif strcmp(rot,'S2C')
    order=[3,2,1];
    res=m*m2;
    
     deg3=-atan2(res(3,2),res(3,3))*180/pi;
   
    m3=rotmat('T2C',deg3);
    
    res=res*m3(:,order);
    
elseif strcmp(rot,'C2S')
    order=[1,3,2];
    res=m*m2;

    deg3=-atan2(res(3,1),res(3,3))*180/pi;
   
    m3=rotmat('T2S',deg3);
    
    res=res*m3(:,order);
    
else
    error('unknown code');
end




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


if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end


res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
