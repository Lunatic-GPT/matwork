
function kern = loadGrid3Kernel(tabsize)	


    kern=zeros(1,tabsize);
    
	for i=2:tabsize-1	
    		kern(i) = kernel(sqrt((i-1)/(tabsize-1)));
    end

    kern(1) = 1.0;

    kern(end) = 0.0;


 function y=kernel(r)
  global    KERNEL_WIDTH
    KERNEL_WIDTH=5;
      OVERSAMPLING_RATIO=1.5;
     
        beta=pi*sqrt((KERNEL_WIDTH/OVERSAMPLING_RATIO*(OVERSAMPLING_RATIO-0.5))^2-0.8);
        y=i0(beta*sqrt(1-r*r))/i0(beta);
        
   function res=i0(x)
            

	ax = abs(x);
	

    if (ax < 3.75) 
    
		y=x/3.75;
        y=y*y;
		res=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492 ...
			   +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	
    else 
    
		y=3.75/ax;
		res=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1 ...
				+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2 ...
				+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1 ...
				+y*0.392377e-2))))))));
    end