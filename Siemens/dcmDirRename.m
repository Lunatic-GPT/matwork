function vargout=dcmDirRename

dir_str=dir('*');
for i=1:length(dir_str)
    
    if strcmp(dir_str(i).name,'.') || strcmp(dir_str(i).name,'..')
        continue;
    end
    
    if exist(dir_str(i).name,'dir')
        
        str=dir([dir_str(i).name,'\*']);
        
        if exist(fullfile(dir_str(i).name,str(end).name),'dir')
         
          str2=dir(  fullfile(dir_str(i).name,str(end).name));
          try
          a=dicominfo(fullfile(dir_str(i).name,str(end).name,str2(end).name));
          
          catch
              continue;
          end
        else
            try
                a=dicominfo(fullfile(dir_str(i).name,str(end).name));
            catch
                % not a dicom dir
                continue;
            end
        end
       
        
        if isfield(a,'SeriesDescription')
            new_dir=strrep(a.SeriesDescription,' ','_');       
            
            n=1;
            while 1
                
                if exist(new_dir,'dir')
                    n=n+1;
                    new_dir=[strrep(a.SeriesDescription,' ','_'),num2str(n)];
                else
                    break;
                end
            end
            
                    
            fprintf('Rename %s to %s;\n',dir_str(i).name, new_dir);
                
            movefile(dir_str(i).name,new_dir);
       
        end
    end
    
end

if nargout>=1
    vargout{1}=dlist;
end

