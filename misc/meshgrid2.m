function [x,y,z]=meshgrid2(x,y,z)
% Meshgrid but use the first dimension as x; The output dim is [length(x),length(y),length(z)]
% In Measgrid, the output dim is [length(y),length(x),length(z)] 
if ~exist('z','var')
    [y,x]=  meshgrid(y,x);
    
    if nargout==1
        x=cat(3,x,y);
    end
else
    [y,x,z]=  meshgrid(y,x,z);
    
    if nargout==1
        x=cat(4,x,y,z);
    end
    
end



