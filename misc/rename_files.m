function rename_files(p,newp,search_sub_dir)
%rename_files(p,newp)

if ~exist('search_sub_dir','var')
    search_sub_dir=false;
end


if strcmp(p,'*')
l=dir2('*');
else
    if ~search_sub_dir
     l=dir(['*',p,'*']);
    else
     l=find_pattern(['*',p,'*']);   
    end
end

for i=1:length(l)
    
    if ~search_sub_dir
    name=l(i).name;
    folder=l(i).folder;
    else
        
        name=filename(l{i});
        folder=fileparts(l{i});
    end
    
    
    if ~strcmp(p,'*')
        pos=strfind(name,p);   
        if isempty(pos)
            continue;
        end
        name2=[name(1:pos-1),newp,name(pos+length(p):end)];
        if exist(name2,'file') || exist(name2,'dir')
            error([name2,' already exists!']);
        end
    else
        pos=strfind(newp,'*');    
        name2=[newp(1:pos-1),name,newp(pos+length(p):end)];
        if exist(name2,'file') || exist(name2,'dir')
            error([name2,' already exists!']);
        end
    end
    disp(['rename ',name,' to ', name2]);
    disp(folder);
    movefile(fullfile(folder,name),fullfile(folder,name2));
    
end

