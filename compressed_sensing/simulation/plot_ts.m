a=BrikLoad('synthesize_kdata_nih_nois0.020_bold0.015+orig');


m=BrikLoad('act_mask2+orig');

bl=[21:40,61:80,101:120,141:160,181:200,221:240];

mn=mean(a(:,:,:,bl),4);
sd=std(a(:,:,:,bl),0,4);

snr=mn./sd;
mean(snr(m(:,:,:,1)>0))
%std(snr(m>0))


%%
ts=a(47,45,4,:);
t=0:2:478;

mn=mean(ts);
figure;plot(t,squeeze(ts)/mn,'-');

y=squeeze(ts)/mn;
ref=load('refg_0_10_30_6cy_TR2.0');

x=[ref(:),ones(length(ref),1)];
b=x\y;

hold on;

plot(t,x*b,'r-');
for i=1:6
    
    x=[(i-1)*80,(i-1)*80,(i-1)*80+20,(i-1)*80+20];
    
    y=[0.95,1.06,1.06,0.95];
    
    
    h=patch(x,y,0.5*ones(1,3));
    set(h,'FaceAlpha',0.5);
end



xlim([0,480]);

ylim([0.95,1.06]);

set(gca,'FontSize',14);

 