function PVSPath(fname,out_prefix)
% m: pvs mask
% prefix: file name of the output file

m=ri_d1(fname);

[c,subind_c]=clusterize2(m>0,2);

tic;
history=sprintf('calculating path for %d PVS\n',length(subind_c));
for i=1:length(subind_c)
    [nvox_path(i),subind{i}]=calcPath(subind_c{i});
end
c_path=c*0;
for i=1:length(subind)
    ind{i}=sub2indb(size(m),subind{i});
    c_path(ind{i}(1:nvox_path(i)))=1;
    c_path(ind{i}(nvox_path(i)+1:end))=2;
end

disp(history);
c=int16(c);
orient=get_orient(fname);
c_path=int16(c_path);
orient.c_path=c_path;
save_pvs(orient,c,ind,subind,nvox_path,out_prefix,history);

