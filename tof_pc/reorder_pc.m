function reorder_pc(dname)

  d=dir(fullfile(dname,'*.dcm'));
for i=1:length(d)
    
    in=dicominfo(fullfile(dname,d(i).name));
    dest=sprintf('%s_%02d.dcm',dname,in.InstanceNumber);
    fprintf('rename %s to %s\n',d(i).name,dest);
    movefile(fullfile(dname,d(i).name),fullfile(dname,dest));
    
end
