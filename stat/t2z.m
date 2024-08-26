function z=t2z(t,dof)
% z=t2z(t,dof)
% converts a t values into z values
%p = tTest(dof,t);
p = tTest(dof,t);
z = sqrt(2)*erfcinv(p)*sign(t);

