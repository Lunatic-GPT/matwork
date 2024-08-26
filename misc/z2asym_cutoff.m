function zasym=z2asym_cutoff(f,z,fw)
% water_ref
% frequency in z spectrum from - to +;
% asym = M(-)-M(+)
% assuming water reference

[f,ind]=sort(f);
z=z(ind);

if fw<=f(1) || fw>=f(end)
    zasym=NaN;
%    error('fw is out of range');
end

if (fw-f(1)) > f(end)-fw
    
 zp=z(f>fw);
 fref=fw-(f(f>fw)-fw);
i1=(f<=fw);
zm=interp1(f(i1),z(i1),fref,'spline');
zasym=zm-zp;

else
    
 zp=z(f<fw);
 fref=fw-(f(f<fw)-fw);
 i1=(f>=fw);
 zm=interp1(f(i1),z(i1),fref,'spline');
 zasym=zm-zp;
end

    
    



