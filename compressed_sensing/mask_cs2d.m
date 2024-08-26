function mask=mask_cs2d(fid_prefix)
%[m2,petable]=mask_cs1d(fid_prefix)
tabfile=readPar(fid_prefix,'petable');

tab=loadtable(tabfile(2:end-1));

cs_pe=readPar(fid_prefix,'cs_pe');
cs_pe2=readPar(fid_prefix,'cs_pe2');
nv=readPar(fid_prefix,'nv');
nv2=readPar(fid_prefix,'nv2');


tab=tab(1:cs_pe*cs_pe2);
mask=zeros(nv,nv2);
mask(tab)=1;



