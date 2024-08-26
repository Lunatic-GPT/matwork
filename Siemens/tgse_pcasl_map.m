function tgse_pcasl_map(dname)

a=ri(dname,1);

dc=mean(a(:,:,2:2:end),3);
dl=mean(a(:,:,3:2:end),3);

a2=(dc-dl)./a(:,:,1);


save(sprintf('perfusion_%s.mat',dname),'a2');

