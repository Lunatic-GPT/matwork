TR=1;
fa=52.7;
T1=2.2;

dm1=(sin(fa*pi/180)-ss_GEEPI(fa,TR,0,T1,0.04));

TR=0.5;
fa=38.8;
dm2=(sin(fa*pi/180)-ss_GEEPI(fa,TR,0,T1,0.04));

dm1/dm2
