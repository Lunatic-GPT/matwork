function J = quadlength_rect(x)
% quadlength -- Find length and dyadic length of square matrix
%  Usage
%    [n,J] = quadlength(x)
%  Inputs
%    x   2-d image; size(n,n), n = 2^J (hopefully)
%  Outputs
%    n   length(x)
%    J   least power of two greater than n
%
%  Side Effects
%    A warning message is issue if n is not a power of 2,
%    or if x is not a square matrix.
%
	s = size(x);
    
    n = s(1);
	k = 1 ; J1 = 0; while k < n , k=2*k; J1 = 1+J1 ; end ;
	if k ~= n
		  disp('Warning in quadlength: n != 2^J');
    end

    
    n = s(2);
	k = 1 ; J2 = 0; while k < n , k=2*k; J2 = 1+J2 ; end ;
	if k ~= n 
		  disp('Warning in quadlength: n != 2^J');
    end
    
    J=[J1,J2];
    
%
% Copyright (c) 1993. David L. Donoho
%     
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
