function stimulate_ijk2pos(fid_prefix,i,j,k,lr_reverse)
%stimulate_ijk2pos(fid_prefix,i,j,k,lr_reverse)
% lr_reverse: whether left and right are switched on the stimulate image.
% default: true

orient=readPar(fid_prefix,'orient');
if ~strcmp('trans90',orient(2:end-1))
    error('only works for axial 90');
end



if exist('lr_reverse','var') && ~lr_reverse
    i=63-i;
end

pro=readPar(fid_prefix,'pro');
ppe=readPar(fid_prefix,'ppe');

lpe=readPar(fid_prefix,'lpe');
lro=readPar(fid_prefix,'lro');

pss=readPar(fid_prefix,'pss');

np=readPar(fid_prefix,'np');
np=np/2;
npe=64;
pos2=-pro+(i-np/2+0.5)*lro/np;
  
pos1=ppe+(j-npe/2+0.5)*lpe/npe;
pos3=pss(k+1);

fprintf('pos1 = %3.2f mm, pos2= %3.2f mm, pos3=%3.2f mm\n',pos1*10,pos2*10,pos3*10);
