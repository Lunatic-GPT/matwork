function clust_name = whereami_roi(mask,tlrc,nvox_min,name_out)
%clust_name = whereami_roi(mask,tlrc,nvox_min)
% mask: mask dataset, same format as in afni
% nvox_min: the minimun number of voxels in a cluster.
if ~exist('nvox_min','var')
    nvox_min = 1;
end
cmd = sprintf('3dclust 0 -%d ''%s'' > 3dclust_tbl.1D',nvox_min,mask);
unix(cmd);

tbl = textread('3dclust_tbl.1D','','commentstyle','shell');
nc = size(tbl,1);
clust_name = cell(1,nc);
fid = fopen(name_out,'w');
for i=1:nc
    abc=orig2tlrc(tbl(i,2),tbl(i,3),tbl(i,4),tlrc);
    cmd = sprintf('whereami %f %f %f 1 -atlas CA_N27_ML',abc(1),abc(2),abc(3));
    [e,m]= unix(cmd);
    o = sep_line(m);
    fprintf(fid,'cluster %d (%3.0f, %3.0f, %3.0f), nvox = %d, mean = %3.2f\n',i,abc(1:3),tbl(i,1),tbl(i,11));
    fprintf(fid,'%s\n',o{:});
    fprintf('cluster %d (%3.0f, %3.0f, %3.0f), nvox = %d, mean = %3.2f\n',i,abc(1:3),tbl(i,1),tbl(i,11));
    fprintf('%s\n',o{:});
    clust_name{i} = o{:};
end
fclose(fid);

function o = sep_line(m)

nl=findstr(m,char(10));

o = {};
for i=1:length(nl)-1
    if ~isempty(findstr(m(nl(i)+1:nl(i+1)-1),'CA_N27_ML'))
        o{end+1} = m(nl(i)+1:nl(i+1)-1);
    end
end
        