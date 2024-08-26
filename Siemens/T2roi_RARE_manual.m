function T2=T2roi_RARE_manual(te,y)
% T2roi_RARE_manual(te,y)

%%
%roi={'roi_csf_T1','roi_gm_T1','roi_wm_T1'};

options=optimset('MaxIter',20,'Display','off');


 [beta,r]=lsqcurvefit(@exp_decay,[max(y),200],te(:),y(:),[],[],options);
 if nargout==0
  figure;
plot(te,y,'o');
 hold on;plot(te,exp_decay(beta,te),'r-');     
 set(gca,'YScale','log');

 title(sprintf('T2 = %s',beta(2)));
 else
 fprintf('T2 = %s\n',beta(2));
     
 end
 T2=beta(2);
