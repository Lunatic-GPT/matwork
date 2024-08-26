function [r,p]=disp_stat(x,dim,title,nanflag)

if exist('title','var')
  fprintf('%s\n',title);
end

if ~exist('nanflag','var')
  nanflag='omitnan';
end

if exist('dim','var') && ~isempty(dim)
 
 fprintf('n=%d;\n',size(x,dim));   
 fprintf('sum = %4.3e; ',sum(x,dim,nanflag));
  disp(' ');
 fprintf('mean = %4.3e; ',mean(x,dim,nanflag));
 disp(' ');
 fprintf('std = %4.3e; ',std(x,[],dim,nanflag));
  disp(' ');
 fprintf('sem = %4.3e; ',std(x,[],dim,nanflag)/sqrt(size(x,dim)));
 disp(' ');

 fprintf('effect size = %4.3e; ',mean(x,dim)./std(x,[],dim));
 disp(' ');
else
 fprintf('n=%d; ',length(x)); 
 fprintf('sum=%4.3e; ',sum(x,nanflag));
 fprintf('mean=%4.3e; ',mean(x,nanflag));
 fprintf('std=%4.3e; ',std(x,nanflag));
 fprintf('sem=%4.3e; ',std(x,nanflag)/sqrt(length(x)));
 fprintf('eff. size=%4.3e; ',mean(x,nanflag)/std(x,nanflag));
 fprintf('range = %4.3e - %4.3e',min(x),max(x));
 disp(' ');    
end


