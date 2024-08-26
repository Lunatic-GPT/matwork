function delete_files_recursive(pattern)
%delete_files(pattern,d_name)



dir_str=dir('*');

mkdir('deleted_files');

for i=1:length(dir_str)
    
    if exist(dir_str(i).name,'dir')
        try
        movefile(fullfile(dir_str(i).name,pattern),'deleted_files\');
        catch
            continue;
        end
    end
end
