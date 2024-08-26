function x=CS_dict(z,m,D,K,thr)

z2=z;

opts.printEvery=0;
opts.maxiter=300;
opts.slowMode=1; 

maxiter=50;
iter=1;


XFM=Wavelet_ISU([size(z2,1),size(z2,2)]);

while iter<=maxiter
 x=fft2c(z2);
 
 for i=1:size(x,1)
     for j=1:size(x,2)
         for k=1:size(x,3)
             b=squeeze(x(i,j,k,:));
             s=OMP(D,b,K,[],opts);
             x(i,j,k,:)=D*s;
             
         end
     end
 end
 
 %%{
 wvl=XFM*x;
 wvl=wthresh(wvl,'h',300);

 x=XFM'*wvl;
 %}
 z2=ifft2c(x);
 z2(m>0)=z(m>0);
 
     
  if exist('x_old','var')    
      dx=x-x_old;
      res=sqrt(sum(abs(dx(:)).^2))/sqrt(sum(abs(x(:)).^2));
      
     fprintf('%d %f\n',iter,res); 
      x_old=x;
      if  res< thr
          break;
      end
      
  else
      x_old=x;
  end
  
  
  iter=iter+1;
end

x=fft2c(z2);


