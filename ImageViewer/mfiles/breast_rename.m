function breast_rename(dlist)
% breast_rename(dlist)



dlist = str2cell(dlist);
c_dur = pwd;
for i=1:length(dlist)
    disp(dlist{i});
    cd(dlist{i});

     
    mkdir('DCE');
    cmd = 'mv *MRDC* DCE';
    unix(cmd);
    
    cd DCE;
    str_dir = dir('*MRDC*');
    cmd = sprintf('re-order %d',length(str_dir));
    unix(cmd);
    
    cd(c_dur);
    cmd = sprintf('chown -R zong:breast %s',dlist{i});
    unix(cmd);
    cmd = sprintf('chmod 770 -R %s',dlist{i});
    unix(cmd);
  
end
    
    



