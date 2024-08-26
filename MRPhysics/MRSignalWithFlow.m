function [ratio,spvs,sart] = MRSignalWithFlow(v,thk,TR,fa,T1art,T1pvs,pvf,sz_init)
%[ratio,spvs,sart] = MRSignalWithFlow(v,thk,TR,fa,T1art,T1pvs,pvf,sz_init)
if ~exist('sz_init','var')
    sz_init=1;
end

spvs=ssFLASH(fa,TR,T1pvs);  % pvs
sart=0;
sz=sz_init;
irep=1;
while irep<=ceil(thk/v/TR)
    sart_tmp=sz*sin(fa/180*pi);
    sz=sz*cos(fa/180*pi);    
    sz=1+(sz-1)*exp(-TR/T1art);   

    if irep==ceil(thk/v/TR)
      sart=sart+sart_tmp*(thk-v*TR*(irep-1))/thk;
    else
      sart=sart+sart_tmp*TR*v/thk;  
    end    
    
    irep=irep+1;
end
    
ratio=(sart*pvf+spvs*(1-pvf))/spvs;

