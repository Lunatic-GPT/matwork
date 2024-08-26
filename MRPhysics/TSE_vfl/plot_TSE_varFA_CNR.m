t=load('TSE_VarFA_CNR_pwr_PVST1T2.mat');


%%
c=squeeze(abs(t.a(1,:,:)-t.a(2,:,:)));

%c(p==0)=0;
iso=(0.05:0.05:0.35)/0.5*100;
c=c/0.5*100;
figure;imshow(c,jet(100));

for is=1:length(iso)
x=[];
y=[];


tmp=find(c(end,1:end-1)<=iso(is)&c(end,2:end)>iso(is));
if ~isempty(tmp)
   jj=tmp(end);
    x(end+1)=jj+(iso(is)-c(end,jj))/abs(c(end,jj+1)-c(end,jj));
   
   y(end+1)=size(c,1);
end

    
for i=1:size(c,2)
    tmp=find(c(1:end-1,i)<=iso(is)&c(2:end,i)>iso(is));
    if ~isempty(tmp)
        x(end+1)=i;
        jj=tmp(end);
        y(end+1)=jj+(iso(is)-c(jj,i))/abs(c(jj+1,i)-c(jj,i));
    end
end


tmp=find(c(1,1:end-1)<=iso(is)&c(1,2:end)>iso(is));
if ~isempty(tmp)
   jj=tmp(end);
    x(end+1)=jj+(iso(is)-c(1,jj))/abs(c(1,jj+1)-c(1,jj));
   
   y(end+1)=1;
end

%hold on;plot(x,y,'w');
end

    
%%
hold on;

iso=[20,40,80,160,300,450,700];

c=t.pwr;

for is=1:length(iso)
x=[];
y=[];


tmp=find(c(end,1:end-1)<=iso(is)&c(end,2:end)>iso(is));
if ~isempty(tmp)
   jj=tmp(end);
    x(end+1)=jj+(iso(is)-c(end,jj))/abs(c(end,jj+1)-c(end,jj));
   
   y(end+1)=size(c,1);
end

    
for i=1:size(c,2)
    tmp=find(c(1:end-1,i)<=iso(is)&c(2:end,i)>iso(is));
    if ~isempty(tmp)
        x(end+1)=i;
        jj=tmp(end);
        y(end+1)=jj+(iso(is)-c(jj,i))/abs(c(jj+1,i)-c(jj,i));
    end
end


tmp=find(c(1,1:end-1)<=iso(is)&c(1,2:end)>iso(is));
if ~isempty(tmp)
   jj=tmp(end);
    x(end+1)=jj+(iso(is)-c(1,jj))/abs(c(1,jj+1)-c(1,jj));
   
   y(end+1)=1;
end

hold on;plot(x,y,'k');
end


%%
set(gca,'Position',[0.2,0.2,0.7,0.7]);
axis on;
ylim([7,51]);
npe=286;
esp=8.34;

ss=0.1:0.1:0.5;
box off;
xtl={};
ind=zeros(1,length(ss));
for i=1:length(ss)
 ind(i)=1+(ss(i)-t.Iss(1))/(t.Iss(end)-t.Iss(1))*(length(t.Iss)-1);
 xtl{i}=sprintf('%2.1f',ss(i));
end
set(gca,'xtick',ind);
set(gca,'XTickLabel',xtl);


ytl={};
TE=200:200:1000;
ind=zeros(1,length(TE));

for i=1:length(TE)
 ind(i)=1+(TE(i)-t.TE(1))/(t.TE(end)-t.TE(1))*(length(t.TE)-1);
 ytl{i}=sprintf('%3.0f',TE(i));
end

set(gca,'ytick',ind);
set(gca,'YDir','Normal');
set(gca,'YTickLabel',ytl);

set(gca,'FontSize',14);
xlabel('Steady Steady Intensity (M_0)');
ylabel('TE (ms)');



set(gca,'Position',[0.1,0.1,0.7,0.7]);