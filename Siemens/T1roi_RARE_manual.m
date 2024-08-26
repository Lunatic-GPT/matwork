function T1roi_RARE_manual(vd,y)
% T1map_sdt(fid_prefix,[t1t2_0])
% fid_prefix.
% t1t2_0: initial values for T1 and T2. default: 2, 0.045
% fit_mean: fit averaged data; default true;
% only tested for single slice 


[vd,ind]=sort(vd);
y=y(ind);


options=optimset('MaxIter',20,'Display','off');


            
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),mean(vd),max(y)],vd(:),y(:),[],[],options);
   %        disp(toc); 
             R1 = 1/beta(2);
            figure;plot(vd,y,'o');
           hold on;plot(vd,exp_decay2(beta,vd),'r-');
           
          title( sprintf('T1 = %f\n',1/R1));
           
           
             
function y=exp_decay2(b,x)
% y=b1*exp(-x/b2)+b3;
% %4.3f*exp(-x/%4.3f)+%4.3f
y=b(1)*exp(-x/b(2))+b(3);


