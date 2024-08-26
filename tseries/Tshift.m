function Tshift(dname)

%
   d=rdSdt(dname);
   pss=readPar(dname,'pss');
  
   [pss_order,ind] = sort(pss);
   
   fd=fft(d,[],4);
   
   sz=size(d);
   for i=1:size(fd,3)
         dt = (ind(i)-1)/length(pss);
       %dt = 1;
         a=exp(-1i*2*pi*dt*((0:sz(4)-1))/sz(4));
         np=floor((length(a)-1)/2);
         a(end-np+1:end)=conj(fliplr(a(2:np+1)));
         
         if mod(length(a),2)==0
             a(length(a)/2+1)=a(length(a)/2+1)*cos(pi*dt);
         end
         a=reshape(a,[1,1,1,sz(4)]);
         fd(:,:,i,:) = fd(:,:,i,:).*repmat(a,[sz(1:2),1,1]);
   end
   
   
   d=ifft(fd,[],4);
   
   writesdt4(d,[dname,'_Tshft']);
   %
   
   
   
   
   
   