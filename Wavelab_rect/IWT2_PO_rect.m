function x = IWT2_PO_rect(wc,L,qmf)
% IWT2_PO -- Inverse 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    x = IWT2_PO(wc,L,qmf)
%  Inputs
%    wc    2-d wavelet transform [n by n array, n dyadic]
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    x     2-d signal reconstructed from wc
%
%  Description
%    If wc is the result of a forward 2d wavelet transform, with
%    wc = FWT2_PO(x,L,qmf), then x = IWT2_PO(wc,L,qmf) reconstructs x
%    exactly if qmf is a nice qmf, e.g. one made by MakeONFilter.
%
%  See Also
%    FWT2_PO, MakeONFilter
%
	Jtmp = quadlength_rect(wc);
	J=min(Jtmp);
    x = wc; 
    
    if Jtmp(1)>=Jtmp(2)
	  nc = [2^(L+1+Jtmp(1)-Jtmp(2)),2^(L+1)];
    else
       nc = [2^(L+1),2^(L+1+Jtmp(2)-Jtmp(1)),]; 
    end
    
	for jscal=L:J-1,
		top = (nc(1)/2+1):nc(1); bot = 1:(nc(1)/2); all = 1:nc(1);
		for iy=1:nc(2)
			x(all,iy) =  UpDyadLo(x(bot,iy)',qmf)'  ...
					   + UpDyadHi(x(top,iy)',qmf)'; 
        end

        top = (nc(2)/2+1):nc(2); bot = 1:(nc(2)/2); all = 1:nc(2);
		for ix=1:nc(1),
			x(ix,all) = UpDyadLo(x(ix,bot),qmf)  ... 
					  + UpDyadHi(x(ix,top),qmf);
		end
		nc = 2*nc;
	end
	
%
% Copyright (c) 1993. David L. Donoho
%     
    
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
