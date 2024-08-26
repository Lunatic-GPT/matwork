function m=matrix_offset_k2tok1(x)
% the matrix to convert from coefficients of spherical harmonics to those
% of polynomials along the six directions.

% 1st dim- x,y,z; 2nd dim - 'xz','yz','z2','xy','x2y2';
x0=x(1);
y0=x(2);
z0=x(3);
m=zeros(3,5);
m(:,1)=[-z0;0;-x0];
m(:,2)=[0;-z0;-y0];
m(:,3)=[x0;y0;-2*z0];
m(:,4)=[-2*y0;-2*x0;0];
m(:,5)=[-2*x0;2*y0;0];


