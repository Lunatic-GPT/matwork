function  res = NUFFT3D(k,w,shift,imSize,nblock,npool)
%res = NUFFT3D(k,w,shift,imSize)
%	non uniform 3D fourier transform 
% 	This is a class which wrappes Jeffery Fesslers
%	NUFFT
%
%	Inputs:
%		k 	- 	n*3 non uniform k-space coordinates (complex values)
%				from -0.5 to 0.5. 
%		w 	-	Density weighting vector
%		shift 	- 	shifts the center of the image [dx,dy,dz]
%		imSize 	- 	size of the image [sx,sy,sz]
%
%	Return:
%		NUFFT object
%
%	Use:
%		FT = NUFFT3D(k,w,shift,imSize);
%		y = FT*x;
%		xx = FT'*y;
% 


if ~exist('nblock','var')
    nblock=1;
end

if ~exist('npool','var')
    npool=1;
end
    res.pools = npool;
    om=k*2*pi;
    
    res.Nd = imSize;
    res.Jd = [5,5,5];
    res.Kd = floor([imSize*1.4]);
    res.n_shift = imSize/2 + shift;
  %  res.st = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');  calculate on
  %  fly to save memory.
    res.nblock=nblock;
    res.om=k*2*pi;
    res.adjoint = 0;
    res.imSize = imSize;
    res.dataSize = [size(k,1),1];
    res.w = sqrt(w);
    
    res = class(res,'NUFFT3D');

