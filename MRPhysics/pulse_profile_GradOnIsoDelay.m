function pulse_profile_GradOnIsoDelay(fname,B1max,pw,offset)
% pulse_profile(fname,B1max,pw,offset)
% fname: bruker pulse name or a vector.
% B1max is in Hz.
% pw is in s.
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

tp=ones(size(rf,1),1)/size(rf,1)*pw;

[tmp,imax]=max(abs(rf));

for j=1:length(offset)
disp(j);

a=dz; 
for i=1:length(rf)
    ind = (i-1)+1;
    
    if i<imax
     H=pulse_H(s,abs(rf(ind)),0,-angle(rf(ind))*180/pi);
    else
     H=pulse_H(s,abs(rf(ind)),offset(j),-angle(rf(ind))*180/pi);
    end
    a=evolve(a,H,tp(i));
 
end 

 my(j)=trace(dy*a);
 mx(j)=trace(dx*a);
 
 mz(j)=trace(dz*a);
 
end

figure;
subplot(2,1,1);plot(offset,sqrt(abs(mx.^2)+abs(my.^2)),'g-');
hold on;
plot(offset,abs(mz),'k-');
legend('M_{perp}','Mz');
title('Magnetization');

xlabel('f (Hz)');

subplot(2,1,2);

plot(offset,phase(mx+1i*my));

xlabel('f (Hz)');


title('Phase');