function make_matlab_shellscript_gpu(fname,func,mem,dur,varargin)
% make_matlab_shellscript(fname,func,mem,dur)
% mem: memory size in Gb.
% duration in hours
% fname: shell script name


fid =fopen(fname,'w');
prefix=strtok2(fname,'.');


    
   
    option=sprintf('srun --ntasks=1 --cpus-per-task=1');
    option=sprintf('%s --time=%d:00:00',option,dur);
    option=sprintf('%s --mem=%dG',option,mem);
    
    option = sprintf('%s --partition=gpu',option);
    option = sprintf('%s --qos gpu_access',option);
    

if ~isempty(varargin)
    if length(varargin)==2
        for i=1:length(varargin{1})
            for j=1:length(varargin{2})
                cmd2=sprintf('%s(%s,%s)',func,num2str(varargin{1}(i)),num2str(varargin{2}(j)));
                fprintf(fid,'%s --job-name %s_%s_%s matbgk "%s" log_%s_%s_%s\n',option,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)),cmd2,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)));
              
            end
        end
    elseif length(varargin)==1
        for i=1:length(varargin{1})
            cmd2=sprintf('%s(%s)',func,num2str(varargin{1}(i)));
            fprintf(fid,'%s --job-name %s_%s matbgk "%s" log_%s_%s\n',option,prefix,num2str(varargin{1}(i)),cmd2,prefix,num2str(varargin{1}(i)));
        end
        
    elseif length(varargin)==3
        for i=1:length(varargin{1})
            for j=1:length(varargin{2})
                for k=1:length(varargin{3})
                    cmd2=sprintf('%s(%s,%s,%s)',func,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)));
                                  fprintf(fid,'%s --job-name %s_%s_%s_%s matbgk "%s" log_%s_%s_%s_%s\n',option,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)),cmd2,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)));
                        
                end
            end
        end
    elseif length(varargin)==4
        for i=1:length(varargin{1})
            for j=1:length(varargin{2})
                for k=1:length(varargin{3})
                    for l=1:length(varargin{4})
                        cmd2=sprintf('%s(%s,%s,%s,%s)',func,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)),num2str(varargin{4}(l)));
                       
                            fprintf(fid,'%s --job-name %s_%s_%s_%s_%s matbgk "%s" log_%s_%s_%s_%s_%s\n',option,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)),num2str(varargin{4}(l)),cmd2,prefix,num2str(varargin{1}(i)),num2str(varargin{2}(j)),num2str(varargin{3}(k)),num2str(varargin{4}(l)));

                    end
                end
            end
        end
        
    end
else

        fprintf(fid,'%s matbgk "%s" log_%s\n',option,func,prefix);       

end

fclose(fid);
