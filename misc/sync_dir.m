function sync_dir(src,dest,flag)
%sync_dir(src,dest,flag)
%flag: 1 - update only + user confirm
%flag: 2 - update only + no prompt
%flag: 99 - update+remove extra + user confirm
if ~exist('flag','var')
    flag=1;
end
print_check(0,0,0);
%fprintf('%10d(%10d) files checked\n',0,0);

s.ntotal=0;
s.file2copy={};
s.file_dest={};
s.dir2rm={};
s.file2rm={};
s.flag=flag;



if iscell(src)
   for i=1:length(src)
    s=compare_dir(src{i},dest{i},s);
   end
else
    s=compare_dir(src,dest,s);
end

if flag<2
   reply = input('update? Y/N[Y]:','s');
else
    reply='Y';
end

if isempty(reply)||strcmp(reply,'Y')
    if flag<2
    reply = input('Are you sure? Y/N[Y]:','s');
    else
       reply='Y'; 
    end
    if isempty(reply)||strcmp(reply,'Y')
        for i=1:length(s.file2copy)
            fprintf('copying %s\n',s.file2copy{i});
            try
            copyfile(s.file2copy{i},s.file_dest{i});
            catch
            fprintf('error copying %s\n',s.file2copy{i});
                
            end
            
            
        end
    end
end

if flag<99
    return;
end
fprintf('extra files in destination\n');
for i=1:length(s.file2rm)
    fprintf('%s\n',s.file2rm{i});
end

for i=1:length(s.dir2rm)
    fprintf('%s\n',s.dir2rm{i});
end

if flag==99
 reply = input('remove extra? Y/N[N]:','s');
else
 reply='Y';
end

if ~isempty(reply)&&(strcmp(reply,'Y')||strcmp(reply,'y'))
if flag==99
    reply = input('Are you sure? Y/N[N]:','s');
else
    reply='Y';
end    
    
    if ~isempty(reply)&&(strcmp(reply,'Y')||strcmp(reply,'y'))
        for i=1:length(s.file2rm)
            delete(s.file2rm{i});
        end
        
        for i=1:length(s.dir2rm)
            rmdir(s.dir2rm{i},'s');
        end
        
    end
    
    
end


function  s=compare_dir(src,dest,s)

a=dir2(src);
if ~exist(dest,'dir')
    mkdir(dest);
end

for i=1:length(a)
    src_fd=fullfile(a(i).folder,a(i).name);
    dest_fd=fullfile(dest,a(i).name);
    
    if isfolder(src_fd)
        s=compare_dir(src_fd,dest_fd,s);
    elseif ~exist(dest_fd,'file')
        s.ntotal=s.ntotal+1;
        s.file2copy{end+1}=src_fd;
        s.file_dest{end+1}=dest_fd;
        nupdate=length(s.file2copy);
        nrm=length(s.file2rm)+length(s.dir2rm);
        print_backspace(48);
              
       % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        fprintf('new file: %s\n',dest_fd);
       % fprintf('%10d(%10d) files checked\n',s.ntotal,nupdate);
        print_check(s.ntotal,nupdate,nrm);
       
    else
        s.ntotal=s.ntotal+1;
        dir_dest=dir(dest_fd);
        
        if a(i).bytes~=dir_dest.bytes || a(i).datenum>dir_dest.datenum
            s.file2copy{end+1}=src_fd;
            s.file_dest{end+1}=dest_fd;
            nupdate=length(s.file2copy);
            nrm=length(s.file2rm)+length(s.dir2rm);
            
       %     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        print_backspace(48);
        
            fprintf('outdated file: %s; size %d vs %d; time %d vs %d\n',...
                dest_fd,a(i).bytes,dir_dest.bytes, a(i).datenum,dir_dest.datenum);
       %     fprintf('%10d(%10d) files checked\n',s.ntotal,nupdate);
        print_check(s.ntotal,nupdate,nrm);
           
            
        end
        
    end
    
    if mod(s.ntotal,1000)==0
        nupdate=length(s.file2copy);
           nrm=length(s.file2rm)+length(s.dir2rm);
         
        %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        print_backspace(48);
        print_check(s.ntotal,nupdate,nrm);
       
       % fprintf('%10d(%10d) files checked\n',s.ntotal,nupdate);
        disp('');
    end
end

if s.flag<99
    return;
end

b=dir2(dest);
for i=1:length(b)
    
            nupdate=length(s.file2copy);
            nrm=length(s.file2rm)+length(s.dir2rm);
            src_fd=fullfile(src,b(i).name);
            dest_fd=fullfile(dest,b(i).name);
    if isfolder(dest_fd) && ~exist(src_fd,'dir')
        %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        print_backspace(48);
        
        fprintf('extra folder: %s\n',dest_fd);
        s.dir2rm{end+1}=dest_fd;
        print_check(s.ntotal,nupdate,nrm);
       % fprintf('%10d(%10d) files checked\n',s.ntotal,nupdate);       
    end
    
    if ~isfolder(dest_fd) && ~exist(src_fd,'file')
        %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); 
        print_backspace(48);
        fprintf('extra file: %s\n',dest_fd);
        s.file2rm{end+1}=dest_fd;
        print_check(s.ntotal,nupdate,nrm);
       % fprintf('%10d(%10d) files checked\n',s.ntotal,nupdate);
       
    end
end


function print_backspace(n)
    for i=1:n
        fprintf('\b');
       % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); 
    end
    
function  print_check(ntotal,nupdate,nrm)
       fprintf('%10d(%10d;%10d) files checked\n',ntotal,nupdate,nrm);
        
        
