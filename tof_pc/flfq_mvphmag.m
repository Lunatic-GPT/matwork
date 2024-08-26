function flfq_mvphmag(maxInd)

if ~exist('maxInd','var')
    maxInd=200;
end
i=1;
while i<=maxInd
    dir_str=dir(sprintf('FL_FQ_*%04d',i));
    
    if ~isempty(dir_str)
        i=i+1;
        
        while 1
            dir_str_p=dir(sprintf('FL_FQ_*P_%04d',i+1));
            dir_str_mag=dir(sprintf('FL_FQ_*MAG_%04d',i));
            
            if ~isempty(dir_str_mag)
                fprintf('move %s to %s\n',dir_str_mag.name,dir_str.name);
                movefile(dir_str_mag.name,dir_str.name);
                i=i+1;
            else
                break
            end
            if ~isempty(dir_str_p)
                
                fprintf('move %s to %s\n',dir_str_p.name,dir_str.name);
                movefile(dir_str_p.name,dir_str.name);
                i=i+1;
            else
                break;
            end
            
        end
    else
        i=i+1;
    end
    
end


