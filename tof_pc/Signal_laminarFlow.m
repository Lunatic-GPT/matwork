 function res=Signal_laminarFlow(r,vmean,va,sa,rad,venc)
       
    vi=v_laminarFlow(r/rad,vmean);
    si=interp1(va,sa,vi);
    res=si.*(exp(1i*vi/venc*pi));
    res(r>rad)=0;



