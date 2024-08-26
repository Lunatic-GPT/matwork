function res = mtimes(a,b)

    
t_dim=ndims(b);      
if a.adjoint
res = 1/sqrt(size(b,t_dim))*fft(b,[],t_dim);

else
 
res = sqrt(size(b,t_dim))*ifft(b,[],t_dim);
end




    
