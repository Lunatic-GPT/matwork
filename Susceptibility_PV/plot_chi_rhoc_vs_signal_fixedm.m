function [pdata,dchi]=plot_chi_rhoc_vs_signal_fixedm(TE,R,B0,theta,dchi_a2)
% vein_QSM_smallAngle_1TE(S,rho,TE,R,B0,theta_B0,dchi_a2); 


dchi=linspace(0.1,1.5,100); 

if ~exist('dchi_a2','var')
    dchi_a2=1.57*(0.0671)^2;
end

if ~exist('B0','var')
   B0=7;
end

if ~exist('theta','var')
    theta = 80;
end

if ~exist('TE','var')
    TE = 14;
end

if ~exist('R','var')
 R = 2; %mm
end

%[ph,maxph]=dchi2phase(0.2,dchi*1e6,B0,TE,0.1,theta,0);
%disp([ph,maxph]);

rhoc=linspace(0,1,100);

s=zeros(length(dchi),length(rhoc));

s_tissue=zeros(length(dchi),length(rhoc));

s_vein=zeros(length(dchi),length(rhoc));

for j=1:length(rhoc)
    for i=1:length(dchi)        
            
        rad=sqrt(dchi_a2/dchi(i));
        [s(i,j),s_tissue(i,j),s_vein(i,j)]=TotalSignal_Vein_Tissue_PerpPlane(rad,dchi(i),B0,TE,theta,rhoc(j),R);
         
    end
end

vox_size=0.4; %mm
vox_snr =10;
nois_vox = vox_size^2/vox_snr;
nois_roi= nois_vox*sqrt(pi*R*R/vox_size^2);

scale_factor=pi*R*R;
sn=s/scale_factor;
s_tissuen=s_tissue/scale_factor;
s_veinn=s_vein/scale_factor;

%%
%% dchi map
cm=jet(100);
pdata=sn;
figure;
    for i=1:length(dchi)
      plot(real(pdata(i,1:end)),imag(pdata(i,1:end)),'-','color',cm(i,:));
      hold on;
    end
colormap(cm);
colorbar(gca);

set(gca,'FontSize',12);
xlabel('Signal Real Part');
ylabel('Signal Imag. Part');


%% rhoc map
% 
% pdata=s_veinn;
% 
% figure;
% for j=1:50%:length(rhoc)
%     
%   %  plot(real(pdata(:,j)),'-','color',cm(j,:));
%        plot(real(pdata(1:51,j)),imag(pdata(1:51,j)),'-','color',cm(j,:));
%        hold on;
% %       plot(real(pdata(21,j)),imag(pdata(21,j)),'or');
% %       plot(real(pdata(47,j)),imag(pdata(47,j)),'or');
% %       plot(real(pdata(48,j)),imag(pdata(48,j)),'or');
% %       plot(real(pdata(74,j)),imag(pdata(74,j)),'or');
% end
% 
% %% s_tissue: map
% 
% figure;plot(dchi*1e6,s_tissuen(:,1),'-');
% set(gca,'FontSize',12);
% xlabel('Susceptibility (ppm)');
% ylabel('Tissue Signal');




