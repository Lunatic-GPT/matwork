function [str_out,nc] = maskcalc(varargin)
%    [str_out,nc] = maskcalc(briks,expr,nvox_min,out_prefix)
% or str_out = maskcalc(briks,expr)
% The first form will generate a mask dataset containing clusters with at
% least nvox_min voxels.  Each cluster is labelled by the size order, starting from 1
% for the largest cluster.
% nc: the number of clusters with nvox>nvox_min in the mask.

if isa(varargin{end-1},'char')
    thr = 1;
    nd = nargin-1;
    clusterized = 0;
    expr = varargin{end};
else
    nd = nargin-3 ;
    clusterized = 1;
    thr = varargin{end-1};
    expr = varargin{end-2};
    nmask = varargin{end};
end

flist = varargin(1:nd);

brik_sel = [];
for i=1:nd
    brik_sel = sprintf('%s -%c %s',brik_sel,'a'+i-1,flist{i});
end
    
dname = sprintf('3dcalc( %s -expr (%s) )',brik_sel,expr);

str_out = dname;

if clusterized

if exist([nmask,'+orig.HEAD'],'file')
    delete([nmask,'+orig.*']);
end
if exist([nmask,'+tlrc.HEAD'],'file')
    delete([nmask,'+tlrc.*']);
end

cmd = sprintf('3dmerge -dxyz=1 -1clust_order 1 %d -prefix %s ''%s''',...
               thr,nmask,dname);
unix(cmd);

if ~isempty(strfind(dname,'+tlrc'))  
   str_out = [nmask,'+tlrc'];
   cmd = ['3dclust 0 -1 ',nmask,'+tlrc.BRIK > 3dclust_tbl.1D'];
   unix(cmd);
else
   str_out = [nmask,'+orig'];   
   cmd = ['3dclust 0 -1 ',nmask,'+orig.BRIK > 3dclust_tbl.1D'];
   unix(cmd);
end


tbl = textread('3dclust_tbl.1D','','commentstyle','shell');
nc = size(tbl,1);
end