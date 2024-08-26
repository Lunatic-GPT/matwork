function M=pulse_profile(fname,B1max,pw,offset)
% pulse_profile(fname,B1max,pw,offset)
% fname: bruker pulse name or a vector.
% B1max is in Hz.
% pw is in s.
% offset in Hz.
    s.nspins=1;
    s.I=1/2;
    s.offset=0;
    
%offset = linspace(-20000,20000,100); 
dy=spin_U(1/2,'y');
dx=spin_U(1/2,'x');
dz=spin_U(1/2,'z');

s.sx=dx;
s.sy=dy;
s.sz=dz;

%[rf,tp]=read_rf_varian('HS-AFP_pi.RF',4e-3,2277);
if isa(fname,'char')
  rf=read_rf_bruker(fname,B1max);
else
  rf=B1max*fname/max(abs(fname));
end

tp=ones(length(rf),1)/length(rf)*pw;

for j=1:length(offset)
disp(j);

a=dz; 
for i=1:length(rf)
    ind = (i-1)+1;
    H=pulse_H(s,abs(rf(ind)),offset(j),-angle(rf(ind))*180/pi);
    a=evolve(a,H,tp(i));
 
end 

 my(j)=trace(dy*a);
 mx(j)=trace(dx*a);
 
 mz(j)=trace(dz*a);
 
end
M=[mx',my',mz'];

figure;
subplot(2,1,1);plot(offset,sqrt(abs(mx.^2)+abs(my.^2)),'g-');
hold on;
plot(offset,mz,'k-');
legend('M_{perp}','Mz');
title('Magnetization');



xlabel('f (Hz)');

subplot(2,1,2);

plot(offset,angle(mx+1i*my));
ph=angle(mx+1i*my);

xlabel('f (Hz)');


title('Phase');