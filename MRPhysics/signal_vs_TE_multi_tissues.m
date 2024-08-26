a=BrikLoad('../../7TData/HVB224/10+orig');
%%
 d=a(:,182,243);
% fd=ifft1c(d,1);
% [d2,f2]=partialFT(fd(81:end)',99);
% 
% 
% figure;plot(d);

gm=[54:66,267:302];
wm=[67:128,159:203,230:266];
csf=[129:158,204:229];

%d(54:302)=1;
d(1:53)=0;
d(303:end)=0;

T2=[87,66,330];
%T2=[160,160,160];

T2d=zeros(size(d));
T2d(gm)=T2(1);
T2d(wm)=T2(2);
T2d(csf)=T2(3);


%TE=[267,370,473,576,3.44*179];
TE=[20,60,120,180,300,500];


tau=3.44;
d3=zeros(length(TE),length(d));
for it=1:length(TE)
    
    nneg=round(TE(it)/tau)-1;
    
    dk=zeros(1,nneg+size(d,1)/2);
    
for i=1:length(dk)
    d2=d.*exp(-i*tau./T2d);
    fd2=ifft1c(d2,1);
    dk(i)=fd2(i+length(d)/2-nneg);
end

d3(it,:)=partialFT(dk,nneg);
    
end


ind={gm,wm,csf};
for i=1:3

    
 opt=statset('FunValCheck','off');
 y=mean(d3(:,ind{i}),2)';

 b=nlinfit(TE,y,@exp_decay,[max(y),100],opt);
 
 fprintf('M0=%f, T2 = %f\n',b(1),b(2));
         
figure;

plot(TE,y,'o');

hold on;
plot(TE,exp_decay(b,TE),'r');
set(gca,'YScale','log');
end
