%% 1D
n=64;
nt=100;
%m2=cstable_1dgauss(n,4,1,nt,1);
%m2=cstable_gems_mTR(n,nt,4,1,6);
m2=cstable_gems_mTR(n,nt,4,1,0);

m=m2(:,8);
m=repmat(m,[1,n]);

figure;
imshow(m,[]);
i1=23;
i2=26;
wl=zeros(n,n);
wl(i1,i2)=1;
xfm=Wavelet_ISU([n,n]);

img=xfm'*wl;


k=fft2c(img);
k2=k;
k2(m==0)=0;
wl2=ifft2c(k2);
wl2=xfm*wl2;
x=-n/2:n/2-1;

[xg,yg]=meshgrid(x,x);
z=abs(wl2)*100;

x2=[x(i2)-1,x(i2),x(i2)+1];
y2=[x(i1)-1,x(i1),x(i1)+1];
[xg2,yg2]=meshgrid(x2,y2);
z2=max(abs(wl2(:)))*zeros(3,3)*100;
z2(2,2)=100*abs(wl2(i1,i2));

figure;
h=mesh(xg,yg,abs(z));
color=gray(200);
set(h,'FaceColor',[1,1,1]);
set(h,'FaceAlpha',1);
colormap((color(50:200,:)));
zlim([0,100]);
hold on;
h2=mesh(xg2,yg2,abs(z2));
set(h2,'FaceColor',[1,0,0]);
set(h2,'FaceAlpha',1);

%%{

%}

xlim([-n/2,n/2-1]);
ylim([-n/2,n/2-1]);
zlim([0,100]);
set(gca,'FontSize',14,'YTick',[-n/2,0,n/2-1],'XTick',[-n/2,0,n/2-1]);







