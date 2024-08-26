function vargout=sss
warning off;
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
        if nargout ==0
            fprintf('%s: ', dir_str(i).name);
        end
        
        if isfield(a,'SeriesDescription')
       %     fprintf('%s;', a.SeriesDescription);
            fprintf(' %s;',a.ProtocolName);
        end
        if nargout ==0
            fprintf('\n');
        end
    end
    
end

if nargout>=1
    vargout{1}=dlist;
end

