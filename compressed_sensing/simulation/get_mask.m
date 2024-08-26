function m=get_mask(pat)

if strcmp(pat,'uniform')

tmp=load('mask_uniform');
m=shiftdim(tmp.m,-1);
m=repmat(m,[64,1,1,1]);

elseif strcmp(pat,'uniform_nc6')
    tmp=load('mask_uniform_nc6');
m=shiftdim(tmp.m,-1);
m=repmat(m,[64,1,1,1]);

elseif strcmp(pat,'gauss')
    
    tmp=load('mask_gauss');
m=permute(tmp.m,[1,2,4,3]);

else
    error('unknown pattern');
end
