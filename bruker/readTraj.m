function a=readTraj(fname)

 fid=fopen([fname,'/traj'],'r','ieee-le');
b=fread(fid,'double');

mtrx=readbPar([fname,'/method'],'PVM_Matrix');
npro=readbPar([fname,'/method'],'NPro');

fclose(fid);

a=reshape(b,[3,mtrx(1)/2,npro]);


