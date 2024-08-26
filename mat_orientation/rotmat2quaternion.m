function [b,c,d,qfac]=rotmat2quaternion(R)

qfac=det(R);
R(:,3)=qfac*R(:,3);

a=0.5*sqrt(1+R(1,1)+R(2,2)+R(3,3));
if a>0
    b=0.25*(R(3,2)-R(2,3))/a;
    c=0.25*(R(1,3)-R(3,1))/a;
    d=0.25*(R(2,1)-R(1,2))/a;
else
    d=sqrt(-(R(1,1)+R(2,2))/2);
    
    if d~=0
    c=R(3,2)/2/d;
    b=R(3,1)/2/d;
    else
     c=sqrt(-(R(1,1)+R(3,3))/2);
     b=sqrt((R(1,1)-R(3,3))/2);
     
     if R(1,2)<0
        b=-b;
     end
        
    end
    
end
