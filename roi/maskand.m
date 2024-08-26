function str_out = maskand(varargin)


brik_sel = [];
expr = [];
for i=1:nargin
    brik_sel = sprintf('%s -%c %s',brik_sel,'a'+i-1,varargin{i});
    if i < nargin
    expr = sprintf('%s%c*',expr,'a'+i-1);
    else
      expr = sprintf('%s%c',expr,'a'+i-1);
    end
end
    
str_out = sprintf('3dcalc( %s -expr ispositive(%s) )',brik_sel,expr);