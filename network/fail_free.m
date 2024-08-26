
function [results,status]=fail_free(cmd,ntrials)
if ~exist('ntrials','var')
    ntrials = Inf;
end

    status=1;
    count = 0;
    results='t';
   while (status~=0 || ~isempty(results)) && count<ntrials  % not work for 3dvolreg and pvsseg.py
       fprintf('...\n');
       cmd2=strrep(cmd,'\','/');
       cmd2=strrep(cmd2,'" "','"\n"');
      
       fprintf(cmd2);
       fprintf('\n');
       [status,results]=system(cmd);
       count = count+1;
      % fprintf('+\n');
    end