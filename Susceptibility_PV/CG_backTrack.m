function x=CG_backTrack(fun,x0,options)

alpha=0.05;
beta=0.6;
step=options.FiniteDifferenceStepSize;
[g,f0]=calcGrad(fun,x0,step);
dx=-g;
x=x0;
iter=0;
fprintf('iter %d: fun = %7.6f\n',iter,f0);

fx=f0;
while iter<options.MaxIterations 

t=1;
  while fun(x+t*dx)>fx+alpha*t*g*dx'
    t=t*beta;
  end
  x=x+t*dx;
  [g2,fx]=calcGrad(fun,x,step);
  
  iter=iter+1;
  fprintf('iter %d: fun = %7.6f; step norm = %7.6f\n',iter,fx,t*sos(dx,2));
  
  if t*sos(dx,2)<options.StepTolerance
      break;
  end
  gamma=sum(g2.^2,2)/sum(g.^2,2);
  dx=-g2+gamma*dx;
  g=g2;
  
  
end



function [G,f0]=calcGrad(fun,x0,step)
f0=fun(x0);

G=0*x0;
for i=1:length(x0)
    dx=0*x0;
    dx(i)=step;
    x=x0+dx;
    G(i)=(fun(x)-f0)/step;    
end

    





