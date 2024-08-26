function res=bp_ana(sys,dia)

hbp=sys>=140 | dia>=90;

prehbp=(sys<140&sys>=120) | (dia>=80&dia<90);

norm=sys<120&dia<80;

res=NaN(size(sys));
res(hbp)=2;
res(prehbp)=1;
res(norm)=0;

