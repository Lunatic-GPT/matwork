
dchi=linspace(0.1,0.45,100)*1e-6; 

m=1.57e-6*(0.0671)^2;

B0=7;
TE = 14;
theta = 80;
rhoc = [0.1,0.3,-1];
R = 1.2; %mm
vox_size=0.4; %mm
vox_snr =10;

%%

s=zeros(length(dchi),length(rhoc));
for j=1:length(rhoc)
    for i=1:length(dchi)        
        
        
        rad=sqrt(m/dchi(i));
        [s(i,j),s_tissue(i,j),s_vein(i,j)]=TotalSignal_Vein_Tissue_PerpPlane(rad,dchi(i),B0,TE,theta,rhoc(j),R);
         
    end
end

nois_vox = vox_size^2/vox_snr;
nois_roi= nois_vox*sqrt(pi*R*R/vox_size^2);
sn=s/nois_roi;

clr=jet(size(sn,2));
figure;
for j=1:size(sn,2)
    hold on;
    for i=1:size(sn,1)
        plot(real(sn(i,j)),imag(sn(i,j)),'.','color',clr(j,:)*i/size(sn,1));
    end
end



figure;
for j=1:size(sn,2)
    hold on;
    plot(dchi,real(sn(:,j)),'-','color',clr(j,:));
    
end



figure;
for j=1:size(sn,2)
    hold on;
    plot(dchi,imag(sn(:,j)),'-','color',clr(j,:));
    
end
legend('1','2','3');
 
 
 