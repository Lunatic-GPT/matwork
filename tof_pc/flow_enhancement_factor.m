function factor= flow_enhancement_factor(v,thk,TR,fa,T1,sz_init)
% factor= flow_enhancement_factor(v,thk,TR,fa,T1,sz_init)
% fa: flig angle in degrees
% calculate enhancement factor (signal ratio with and without flow) due to flow
% flow enhancement factor for a square pulse

if ~exist('sz_init','var')
    sz_init=1;
end

%T1art=0.2;
%T1tissue=1.4;


ss=ssFLASH(fa,TR,T1);  % WM
sf=0;
sz=sz_init;


if v==0
    factor = 1;
    return;
end

sart_tmp=zeros(1,ceil(thk/v/TR));

for irep=1:ceil(thk/v/TR)
    sart_tmp(irep)=sz*sin(fa/180*pi);
    sz=sz*cos(fa/180*pi);    
    sz=1+(sz-1)*exp(-TR/T1);   
end



for irep=1:ceil(thk/v/TR)
    
    if irep==ceil(thk/v/TR)
      sf=sf+sart_tmp(irep)*(thk-v*TR*(irep-1))/thk;
    else
      sf=sf+sart_tmp(irep)*TR*v/thk;  
    end
    
end
     
factor=sf/ss;

