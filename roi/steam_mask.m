function steam_mask(steam_fid_prefix,parent_grid,out_prefix)
% steam_mask(steam_fid_prefix,parent_grid,out_prefix)
% out_prefix: default 'mask_steam'

if ~exist('out_prefix','var')
    out_prefix='mask_steam';
end

name = steam_fid_prefix;
vox1=readPar(name,'vox1');
vox2=readPar(name,'vox2');
vox3=readPar(name,'vox3');
pos1=readPar(name,'pos1');
pos2=readPar(name,'pos2');
pos3=readPar(name,'pos3');
voxpos2mask([vox2,vox1,vox3],[pos2,pos1,-pos3]*10,parent_grid,out_prefix);
