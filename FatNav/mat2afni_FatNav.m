function mat2afni_FatNav(prefix,lstart,center)
d=single([]);

for i=1:length(lstart)
fname=sprintf('%s_%d_*.mat',prefix,lstart(i));
   disp(fname); 
fname=dir(fname);
 
    fname=fname.name;
    dtmp=load(fname);
    d=cat(4,d,dtmp.d);

end

%% sagital navigator echo
in.ORIENT_SPECIFIC = [5,2,1];  

sz=size(d);

in.DELTA = [-2,-2,-2];
in.ORIGIN= -in.DELTA.*sz(1:3)/2+center;

WriteBrikEZ(d,in,prefix,prefix);
