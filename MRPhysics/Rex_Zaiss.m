function rex=Rex_Zaiss(R2b,fb,kb,b1,Dba,delta)
% res=Rex_Zaiss(R2b,fb,kb,b1,Dab,delta)
% Dba is the chemical shift between pool b and a in Hz.
% delta is the frequency difference between rf and pool a in Hz.


w1=b1*2*pi;
Dba=Dba*2*pi;
delta=delta*2*pi;

Db=delta-Dba;
Gamma=2*sqrt((kb+R2b)/kb*w1^2+(kb+R2b)^2);

th=atan(w1/delta);
rexmax=fb*kb*sin(th)^2*(Dba^2+R2b/kb*(w1^2+delta^2)+R2b*(kb+R2b))/Gamma^2*4;
rex=rexmax*Gamma^2/4/(Gamma^2/4+Db^2);
