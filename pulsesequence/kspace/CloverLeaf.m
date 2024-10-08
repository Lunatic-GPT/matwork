function res = CloverLeaf(circ_scale)

%CloverLeaf(200,200,90,0.1);

GRAD_RASTER_TIME = 1e-5;
nramp2nrad=1;
ncirc2nrad=0.45;

BwPerPixel = 400; %Hz
FOV = 220; % mm
VoxSize = 0.5; % mm
nvox=round(FOV/VoxSize);
ind = zeros(6,2);  % start and end indices for the k-space tracks; ky; x-y arc; kx; x-z arc; kz; y-z arc 

BwPerPixel=BandWidthPerPixel(BwPerPixel,nvox,1);

dt=1/(BwPerPixel*nvox);  %s
np_GRT = GRAD_RASTER_TIME/dt;

if np_GRT>=1
    nrad=round(nvox/np_GRT/2);
    nvox=nrad*np_GRT*2;
else
    nvox=round(nvox/2)*2; 
    nrad=nvox/np_GRT/2;     
end

FOV = nvox*VoxSize;

kmax=pi/VoxSize; % 1/mm

narc=round(nrad*pi/2);
nramp=round(nrad*nramp2nrad);
ncirc=round(nrad*ncirc2nrad);

lstart= -1/3;  % start of k-space position for constant gradient

nneg=round(-nrad*lstart);


if np_GRT<1
    narc = round(narc*np_GRT)/np_GRT;
    ncirc = round(ncirc*np_GRT)/np_GRT;
    nneg = round(nneg*np_GRT)/np_GRT;
    
end

phi=linspace(0,pi/2,narc+1);
xarc=cos(phi(2:end))';
yarc=sin(phi(2:end))';

[ycirc,xcirc]=traj_return(ncirc+1,1/nrad/circ_scale);

xcirc=xcirc(2:end)*circ_scale;
ycirc=ycirc(2:end)*circ_scale;
xcirc=xcirc(:);
ycirc=ycirc(:);

xrad=linspace(0,1,nrad+1);
xrad=xrad(2:end)';

famp=@(x) 2*lstart/(cos(nramp*2*pi/x)-1);

fitfunc=@(x) famp(x)/2*sin(nramp*2*pi/x)*2*pi/x*nrad+1; 

nT  = fsolve(fitfunc,nramp*2);
amp = famp(nT);

ar = [zeros(nramp,1),(cos((1:nramp)'*2*pi/nT)-1)*amp/2,zeros(nramp,1)];


aneg=linspace(lstart,0,nneg+1);

aneg=aneg(2:end)';

aneg=[0*aneg,aneg,0*aneg];
a = [0*xrad,xrad,0*xrad];

a2 = [-ycirc,xcirc+1,0*ycirc];

a3=[yarc,xarc,0*xarc(:)];


a4=[ycirc+1,-xcirc,0*ycirc];


a5=[1-xrad,0*xrad(:),0*xrad(:)];
a6=[-xrad,0*xrad(:),0*xrad(:)];
%a7=[

b=cat(1,ar,aneg,a,a2,a3,a4,a5,a6);  %xy plane


ind_1c=size(ar,1)+size(aneg,1)+1;
ind_rostart=size(ar,1)+1;
ind(1,1)=1;
ind(1,2)=ind(1,1)+size(a,1)-1;
ind(2,1)=ind(1,2)+size(a2,1)+1;
ind(2,2)=ind(2,1)+size(a3,1)-1;
ind(3,1)=ind(2,2)+size(a4,1)+1;
ind(3,2)=ind(3,1)+size(a5,1)+size(a6,1)-1;


figure;plot(b(:,1),b(:,2));
xlim(min_max(b(:,1)));
ylim(min_max(b(:,2)));

a7 = [-xcirc-1,0*ycirc,ycirc];
a8= [-xarc,0*yarc,-yarc];
a9=[xcirc,0*ycirc,-ycirc-1];
a10=[0*xrad,0*xrad,xrad-1];
a11=[0*xrad,0*xrad,xrad];

b2=cat(1,a7,a8,a9,a10,a11); % xz plane;

ind(4,1)=ind(3,2)+size(a7,1)+1;
ind(4,2)=ind(4,1)+size(a8,1)-1;
ind(5,1)=ind(4,2)+size(a9,1)+1;
ind(5,2)=ind(5,1)+size(a10,1)+size(a11,1)-1;



% figure;plot(b2(:,1),b2(:,3));
% xlim([-1.2,1.2]);
% ylim([-1.2,1.2]);


a12 = [0*xcirc,ycirc,xcirc+1];
a13=[0*xarc,-yarc,xarc];
a14=[0*xcirc,-ycirc-1,-xcirc];

ind(6,1)=ind(5,2)+size(a12,1)+1;
ind(6,2)=ind(6,1)+size(a13,1)-1;

tt = round(nrad*1.2);
amp=1/nrad;
tf=round((1-amp/2*tt)*2/amp);

a15=amp*(1:tf)'-1;
a15=[0*a15,a15,0*a15];
tr=tt-tf;
a16 = a15(end,2)+amp/tr*(tr*(1:tr)-(1:tr).^2/2)';
a16=[0*a16,a16,0*a16];


b3=cat(1,a12,a13,a14,a15,a16); % xz plane;
% figure;plot(b3(:,2),b3(:,3));
% 
% xlim([-1.2,1.2]);
% ylim([-1.2,1.2]);

res=cat(1,b,b2,b3);

gamma=42580*2*pi; %Hz/mT

k=res*kmax; % 1/mm

Gconst = BwPerPixel*nvox/gamma*2*pi/FOV; % mT/mm

G = diff(res,1,1);

G1 = G*Gconst/G(ind_1c,2)*1000;  % mT/m
G1(end,:)=0;
% or G can be calculated from the difference of k
% k*x=G*x*gamma*dt; 

%G2 = diff(k,1,1)/nvox*nrad*2/gamma/dt*1000;  % mT/m

figure;

for i=1:3
    subplot(4,1,i);
    plot(k(:,i));
end
subplot(4,1,4);
plot(sos(k,2));

figure;

for i=1:3
    subplot(4,1,i);
    plot(G1(:,i));
end
subplot(4,1,4);
plot(sos(G1,2));

figure;

dt2=dt*nvox/nrad/2;  % gradient wave form time unit 

for i=1:3
    subplot(3,1,i);
    plot(diff(G1(:,i))/GRAD_RASTER_TIME*0.001);
if i==1
    title('Slew rate');
end
ylabel('mT/m/ms');
end



%fprintf('Total Time = %f\n',size(G2,1)*GRAD_RASTER_TIME);

disp('Gradient Size');
disp(size(G1));
fid=fopen('grad3.txt','w');
fprintf(fid,'%d,%d,%d,%d,%d\n',size(G1,1),ind_1c,ind_rostart,nrad,dt*1e9);

for i=1:6
fprintf(fid,'%d,%d\n',ind(i,1),ind(i,2));
end

fprintf(fid,'%f\n',max(abs(G1(:))));

G1=G1/max(abs(G1(:)));

for i=1:3
    
   % fprintf(fid,'{');
 
    for j=1:size(G1,1)
       fprintf(fid,'%4.3f,',G1(j,i));
       if mod(j,50)==0
           fprintf(fid,'\n');
       end
    end
    
end


fclose(fid);

    




