function not_found=wait_for_files(prefix,ind1,ind2,nsecond,do_wait)
found =false;

if ~exist('do_wait','var')
    do_wait=true;
end

if ~exist('nsecond','var')
    nsecond = 60;
end

while ~found
    pause(nsecond);
    found=true;
    
    if isempty(ind1)
         name=prefix;
            
            if ~exist(name,'file')
                fprintf('%s not found\n',name);
                found = false;
            else
                try 
                    whos('-file',name);  % file not ready, although exist
                catch
                    found = false;
                end
                
            end
    else
        for i=1:length(ind1)
            
            name=sprintf('%s_%d_%d.mat',prefix,ind1(i),ind2(i));
            
            if ~exist(name,'file')
                fprintf('%s not found\n',name);
                found = false;
                break;
            else
                try 
                    whos('-file',name);  % file not ready, although exist
                catch
                    found = false;
                    break;
                end
                
            end
            
        end
    end
    
    if ~do_wait
        break;
    end
end

not_found=~found;
