function calcShortestFlowComp%(lRD,lFTaC,dAmp,lRU,venc,Gmax,F1,G1,F2,G2)

%syms RD RU F1 G1 F2 G2 Gs Fa;
% when not inverted
Gs=12.033079e-1;
Fa=1015e-6;
RDs=120e-6;
Gmax=3.4;
%% when inverted and venc=-1
RD = 330e-6;
RU = 330e-6;
venc=4;
F1=2010e-6;
G1=3.4678;

F2=2380e-6;
G2=-3.4716696;
m1_venc=venc2M1(venc);

m1=calc_m1(RD,RU,F1,G1,F2,G2,Gs,Fa,RDs);
disp([m1,m1_venc]*1000);

venc=4;
if venc<0
    G1=Gmax;
    G2=-Gmax;
else
    G1=-Gmax;
    G2=Gmax;
end

a=- 12*G1*G2 + 12*G2*G2;
b= - 24*G1*G2*RD - 24*G1*G2*RU + 12*G2*G2*RD + 12*G2*G2*RU;
c=  12*Fa*Fa*G1*Gs - 12*Fa*Fa*Gs*Gs + 24*Fa*G1*Gs*RDs + 12*Fa*G1*Gs*RU - 12*Fa*Gs^2*RDs + 11*G1*G1*RD*RD + G1*G1*RU*RU + 2*G1*G2*RD*RD - 18*G1*G2*RD*RU - 8*G1*G2*RU*RU + 20*G1*Gs*RDs*RDs + 6*G1*Gs*RDs*RU + 3*G2*G2*RD*RD + 6*G2*G2*RD*RU + 3*G2*G2*RU*RU - 3*Gs*Gs*RDs*RDs;

m1_venc=venc2M1(venc);
a=-a/(24*G1);
b=-b/(24*G1);
c=-c/(24*G1)-m1_venc;
F2=(-b+sqrt(b*b-4*a*c))/2/a;
if F2<0
    F2=(-b-sqrt(b*b-4*a*c))/2/a;
end

F1=-(F2*G2 + Fa*Gs + (G1*RD)/2 + (G2*RD)/2 + (G1*RU)/2 + (G2*RU)/2 + (Gs*RDs)/2)/G1;

fprintf('%6.2f ',[F1,F2]*1e6);
fprintf('\n');
M0_1=G1*(F1+RD);
M0_2=G2*(F2+RD);

F1=ceil(F1*1e5)/1e5;
F2=ceil(F2*1e5)/1e5;

G1=M0_1/(F1+RD);
G2=M0_2/(F2+RD);


fprintf('%6.2d ',[F1,F2]*1e6);
fprintf('\n');
m0=calc_m0(RD,RU,F1,G1,F2,G2,Gs,Fa,RDs);
fprintf('m0 = %f\n',m0);
m1=calc_m1(RD,RU,F1,G1,F2,G2,Gs,Fa,RDs);
disp([m1,m1_venc]*1000);

function m0=calc_m0(RD,RU,F1,G1,F2,G2,Gs,Fa,RDs)
m0=Gs*Fa;
m0=m0+M0Ramp(RDs,Gs);
m0=m0+M0Trap(RU,RD,F1,G1);
m0=m0+M0Trap(RU,RD,F2,G2);
%m0= simplify(m0);

function m0=M0Trap(RU,RD,F,G)

m0=M0Ramp(RU,G);
m0=m0+G*F;
m0=m0+M0Ramp(RD,G);

function m1=calc_m1(RD,RU,F1,G1,F2,G2,Gs,Fa,RDs)

m1=M1Flat(0,Fa,Gs);
t0=Fa;
m1=m1+M1Rampdown(t0,RDs,Gs);
t0=t0+RDs;

m1=m1+M1Trap(t0,RU,RD,F1,G1);
t0=t0+RU+RD+F1;
m1=m1+M1Trap(t0,RU,RD,F2,G2);

%m1= simplify(m1);

function m1=M1Trap(t0,RU,RD,F,G)

m1=M1Rampup(t0,RU,G);
t0=t0+RU;
m1=m1+M1Flat(t0,F,G);
t0=t0+F;
m1=m1+M1Rampdown(t0,RD,G);


function m1=M1Flat(start, D,G)

m1=G*(start+D)^2/2-G*start^2/2;


function m1=M1Rampup(t0, r,G) % Gradient is zero at t0 and G at t0+r
% t0 is the start time;

m1=G/r*( (t0+r)^3/3-t0*(t0+r)^2/2  )-G/r*( t0^3/3-t0^3/2  );

function m1=M1Rampdown(t0, r,G) % Gradient is G at t0, and 0 at t0+r
% t0 is the start time;
S=-G/r;

m1=S*( t0^3/3-t0^3/2)  - S*( (t0-r)^3/3-t0*(t0-r)^2/2 );

function m0=M0Ramp(r,G)
m0=G*r/2;

