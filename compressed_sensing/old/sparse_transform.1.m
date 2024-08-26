function fm=sparse_transform(m,us,dagger)


if ndim(m)<3
    if ~dagger
     fm=fft2(m);
     fm(~us)=0;
    else
     fm = fft2(m.');
     fm = fm';
     fm(~us)=0;
    end
    
else
    
    error('not supported yet');
end
    
    

