function resample_roi(m,n1,n2)

n=lcm(n1,n2);

prefix=strtok(m,'.');

try
m=ri(m,'','','d');
catch
m=ri(m);
    
end
m=repmat2(m,n/n1);

sk=n/n2;
m=m(1:sk:end,1:sk:end,:);


save(sprintf('%s_%dTo%d.mat',prefix,n1,n2),'m');



