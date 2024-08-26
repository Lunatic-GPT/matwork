ts=[ones(1,10),zeros(1,10)];
ts=repmat(ts,[1,4]);
ref=ts;
p=zeros(1,10000);
for i=1:10000
%    ts=randn_correlated(size(ref));
    ts=randn_white(size(ref)); 
    [cc,ptmp]=corrcoef(ts,ref);
    
    p(i)=ptmp(1,2);
    disp(i);
end

figure;plot(ts);
ylim([-0.5,1.2]);