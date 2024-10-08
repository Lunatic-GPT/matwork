function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint
    res = reshape(b,a.dataSize);
  %  res = zpad(bb.*a.mask,a.imSize(1),a.imSize(2));
    res=res.*a.mask;
    
      res = ifft2c(res);
      res=ifft1c(res,3);  
      res = res.*conj(a.ph);
    switch a.mode
    	case 0
		res = real(res);
   	case 1
		res = real(res);
    end



else
    bb = reshape(b,a.imSize);
    
    switch a.mode
    	case 0
		bb = real(bb);
   	case 1
		bb = real(bb);
    end
    
    bb = bb.*a.ph; % phase correct
    res = fft2c(bb);
    res=fft1c(bb,3);
    %res = crop(res,a.dataSize(1),a.dataSize(2));
    res = res.*a.mask;
end

if size(b,2) == 1
    res = res(:);
end



    
