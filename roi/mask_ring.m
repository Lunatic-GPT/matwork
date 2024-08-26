function m=mask_ring(dim,rad,center,include_equal)
% m=mask_circle(dim,rad,center,include_equal)
% rad: 1*2; the inner and outer diameters

 r1=mask_circle(dim,rad(1),center,include_equal);
    r2=mask_circle(dim,rad(2),center,include_equal);
    
    m=r2>0&r1==0;
    
