function T1map_RARE_Siemens_dcm(d)
% T1map_sdt(fid_prefix,[t1t2_0])
% fid_prefix.
% t1t2_0: initial values for T1 and T2. default: 2, 0.045
% fit_mean: fit averaged data; default true;
% only tested for single slice

img2=ri(d,1);
%vd=readdPar(d,'RepetitionTime',1);
extp(d);
longTR = readsPar([d,'.pro'],'adFree[4]');
    shortTR = readsPar([d,'.pro'],'adFree[5]');
    nTR = readsPar([d,'.pro'],'alFree[5]');
    
    TI=logspace(log10(shortTR),log10(longTR),nTR);
    vd=TI/1000;
    
options=optimset('MaxIter',20,'Display','off');

tmp=img2(:,:,:,end);
m=tmp>0.02*max(tmp(:));
R1=0*img2(:,:,:,1);
res1=R1;
sz=size(img2);
for i=1:sz(1)
    fprintf('%d/%d\n',i,sz(1));
    for j=1:sz(2)
        for k=1:sz(3)
          if m(i,j,k)==0
              continue;
          end
            
            
            y = squeeze(img2(i,j,k,:));
        
            [beta,r]=lsqcurvefit(@exp_decay2,[-max(y),mean(vd),max(y)],vd(:),y(:),[],[],options);
   %        disp(toc); 
            if ~any(isnan(beta)) 
             ss = sum((y-mean(y)).^2);   
             R1(i,j,k) = 1/beta(2);
             res1(i,j,k) = 1-sum(r.^2)/ss;
            end
    
        end
           
    end
end

  save([d,'_R1'],'R1','res1');


function y=exp_decay2(b,x)
% y=b1*exp(-x/b2)+b3;
% %4.3f*exp(-x/%4.3f)+%4.3f
y=b(1)*exp(-x/b(2))+b(3);


