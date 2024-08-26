function mean_wght_ntrls(varargin)

ntrl = zeros(nargin-1,0);

for i=1:nargin-1
    
 ntrl(i) = readPar(varargin{i},'arraydim');

end

