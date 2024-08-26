function compare_pair(varargin)
%compare_pair(g1,g2,test)
%or compare_pair(g,test): g dimension:n*2

if nargin==2
    g1=varargin{1}(:,1);
    g2=varargin{1}(:,2);
    test=varargin{2};
    
elseif nargin==3
    g1=varargin{1};
    g2=varargin{2};
    test=varargin{3};
else
    error('number of argumnets wrong');
end

pct=mean(g1./g2)*100-100;
epct=std(g1./g2)*100;

if strcmp(test,'signrank')
  pval=signrank(g1-g2);
  pval_pct=signrank(g1./g2-1);
else
    pval=ttest_1vec(g1-g2);
    
    pval_pct=ttest_1vec(g1./g2-1);
end

   fprintf('%4.3f (%4.3f) vs %4.3f (%4.3f); p = %4.3g;  N = %d\n',mean(g1),std(g1),...
        mean(g2),std(g2),pval,length(g1)); 
    
   fprintf('pct diff = %4.3f (%4.3f); p = %4.3f\n',pct,epct,pval_pct); 
    