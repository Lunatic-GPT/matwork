function [y,x]=traj_return(ncirc,init_stepx)



xx = linspace(0,pi/2,ncirc); % desired sampling


if ncirc<=100

x1=xx(2);
x2=xx(end-1);

else

x1=xx(3);
x2=xx(end-2);
init_stepx=init_stepx*2;    
    
end
x=[0,x1,pi/4,x2,pi/2];
r=[0,init_stepx/cos(x1),1,init_stepx/cos(x1),0];




k = 4; 

sp = spapi( optknt(x,k), x, r );


yy = fnval(xx,sp);

y=yy.*sin(xx);
x=yy.*cos(xx);
if nargout==0
figure;plot (x,y); %
xlim([0,0.8]);
ylim([0,0.8]);

end