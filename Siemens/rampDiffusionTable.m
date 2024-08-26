function rampDiffusionTable(table, rampTR,pwr)

a=read_diffvectors(table);

scl=(0:size(a,1)-1)/rampTR;

scl(scl>1)=1;

scl=scl.^pwr;

b=a.*repmat(scl(:),[1,3]);

prefix=strtok(table,'.');
prefix=sprintf('%s_RampPwr%2.1f',prefix,pwr);
write_DiffusionVectors(b,prefix);







