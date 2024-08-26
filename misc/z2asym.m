function [zasym,fz]=z2asym(f,z,fw)
% water_ref
% frequency in z spectrum from - to +;
% asym = M(-)-M(+)
% assuming water reference

[f,ind]=sort(f);
z=z(ind);


zp=z(f>fw);
zm=zeros(size(zp));

fref=fw-(f(f>fw)-fw);

nf=length(find(fref>=f(1)));
i1=(f<=fw);
zm(1:nf)=interp1(f(i1),z(i1),fref(fref>=f(1)),'linear');
zm(nf+1:end)=z(1);

zasym=zm-zp;
fz=f(f>fw)-fw;





