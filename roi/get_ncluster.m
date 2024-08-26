function nc = get_ncluster(mask)
cmd = ['3dclust 0 -1 ',mask,' > 3dclust_tbl.1D'];
unix(cmd);
tbl = textread('3dclust_tbl.1D','','commentstyle','shell');
nc = size(tbl,1);