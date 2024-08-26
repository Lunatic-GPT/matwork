function fm=sparse_FT(m,us,dagger)


if ndims(m)<3
    if ~dagger
     fm=fft2(m);
     fm(~us)=0;
    else
     m(~us)=0;   
     fm = fft2(conj(m));
     fm = conj(fm);
    end
    
else
    
    error('not supported yet');
end
    
    

