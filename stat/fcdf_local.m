function d=fcdf_local(f,dof1,dof2)
% this definition conflict with a standard matlab function 6/30/10   
x = dof1*f/(dof1*f+dof2);
d = betainc(x,dof1/2,dof2/2)/beta(dof1/2,dof2/2);
