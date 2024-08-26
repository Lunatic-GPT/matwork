function res=mat_SLICE_POS(orient,deg)
% reference to Cordinate_system.txt for convention
% output  y z are opposite of dicom convention
%res: first columne pe; 2nd ro; 3rd 3D
m=eye(3);
deg=(deg-90)/180*pi;

if strcmp(orient,'C')
    norm=SliceOrient2Norm('C2T',0);  %dicom convention

      m2d=[cos(deg),sin(deg);-sin(deg),cos(deg)];
      m(1:2,1:2)=-m2d;
   %   norm(2:3)=-norm(2:3);
  %  rotmat([1,3],[1,3])=m2d;
      
elseif strcmp(orient,'T')
     norm=SliceOrient2Norm('T2C',0);  %dicom convention
     m2d=[cos(deg),sin(deg);-sin(deg),cos(deg)];
     m(1:2,1:2)=m2d;
    
elseif strcmp(orient,'S')
      norm=SliceOrient2Norm('S2T',0);  %dicom convention
      m2d=[cos(deg),sin(deg);-sin(deg),cos(deg)];
      m(1:2,1:2)=m2d;
    
   %   rotmat([3,2],[3,2])=m2d;
else
    
    norm=orient;
    ori=Norm2SliceOrient(norm);
    if strcmp(ori(1),'S')
      m2d=[cos(deg),sin(deg);-sin(deg),cos(deg)];
      m(1:2,1:2)=m2d;
    end
    res=Norm2InplaneDirection(norm);
    
   
end

 [theta,phi]=unitVec2thetaPhi(norm(:));
 res=transform_matrix_rotation(theta,phi);
% disp([det(res),det(m)]);
 res=res*m;
  res(2:3,:)=-res(2:3,:);
 
 
 
 
  
 