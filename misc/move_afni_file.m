function move_afni_file(f1,f2)

f2_prefix = strtok(f2,'+');
f1_prefix = strtok(f1,'+');

      if ~exist([f1_prefix,'+orig.HEAD'],'file')
          error([f1_prefix,'+orig.HEAD'], ' does not exist');
      end
           
      if exist([f2_prefix,'+orig.HEAD'],'file')
       delete(sprintf('%s+orig*',f2_prefix));
      end
      
       cmd = sprintf('3dcopy %s+orig %s', f1_prefix,f2_prefix);
       unix(cmd);
       
       delete([f1_prefix,'+orig*']);
       