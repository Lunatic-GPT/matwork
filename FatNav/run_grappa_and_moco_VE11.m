
function run_grappa_and_moco_VE11(fname)
% sub: subject number (int)
% fname: mat file name
% assumed data already in the current directory and pro file saved as
% prefix/prefix.pro
% motion file under ../../motionFile/
run('~/matlab/fessler/setup.m');

mid=strtok_no(fname,'_',2);

kdata=recon_tse_vfl_Siemens_VE11(fname,['Grappa_',mid]);
%kdata=zeros(2,300,208,2);
dfile=sprintf('motion_%s.1D',mid);
prefix_out=sprintf('Moco_%s',mid);
moco_tse_vfl_afterGrappa_VE11(kdata,fname,dfile,8,prefix_out);