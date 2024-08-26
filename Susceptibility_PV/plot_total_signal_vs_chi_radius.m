
dchi=linspace(0,1,100)*1e-6; 



%m=0.4e-6*0.1*0.1;

B0=7;
TE = 14;
theta = 90;
rhoc = [0.1,0.3,0.5,0.7,0.9];
R = 4; %mm
vox_size=0.4; %mm
vox_snr =10;

rad=linspace(0.01,vox_size,50);
s=zeros(length(dchi),length(rad));
for k=1:length(rhoc)
for j=1:length(rad)
    for i=1:length(dchi)
        [s(i,j),s_tissue(i,j),s_vein(i,j)]=TotalSignal_Vein_Tissue_PerpPlane(rad(j),dchi(i),B0,TE,theta,rhoc(k),R);
        
        
        
        %[s(i,j)=TotalSignal_Tissue_PerpPlane(sqrt(m/dchi(i)),dchi(i),7,10,90,1,0.3);
    end
end

nois_vox = vox_size^2/vox_snr;
        
nois_roi= nois_vox*sqrt(pi*R*R/vox_size^2);
        
sn=s/nois_roi;

clr=jet(size(sn,2));
figure;
for j=1:size(sn,2)
 hold on;
 plot(real(sn(:,j)),imag(sn(:,j)),'color',clr(j,:));
end

title(sprintf('Vein Signal = %2.1f\n',rhoc(k)));
end




 
 
 