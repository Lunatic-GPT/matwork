function R=quaternion2rotmat(b,c,d)

a=sqrt(1-b*b-c*c-d*d);

 R(1,:)= [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ];
 R(2,:) = [ 2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ];
 R(3,:)= [ 2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ];
         
