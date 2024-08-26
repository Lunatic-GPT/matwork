function [dchi_a2,resid,cc,interc]=magMoment_outSidePattern(ph,TE,B0,th_B0,mask,center,voxSize,plot_results)

% B0 projection is along the second dimension
% center: indices for the center 
% dchi_a2 unit: ppm*mm2
% th_B0 and ph in degree
% TE in ms;
% B0 in tesla

if ~exist('plot_results','var')
    plot_results=false;
end
if numel(voxSize)==1
    voxSize=voxSize*ones(1,length(center));
end

TE = TE/1000;
gamma=42.58e6*2*pi;

pos=find(mask>0);
pos=ind2subb(size(mask),pos);
ndim=ndims(mask);
pos=pos-repmat(center,size(pos,1),1);
pos=pos.*repmat(voxSize(1:ndim),size(pos,1),1);

phi=atan2(pos(:,1),pos(:,2));

x=cos(2*phi).*TE*sin(th_B0/180*pi).^2*B0*gamma/2./sum(pos(:,1:2).^2,2)/1e6;


y=ph(mask>0)*pi/180; 


xx=[x(:),ones(length(x),1)];

b=xx\y;
  
dchi_a2=b(1);
interc=b(2);

resid=sos(xx*b-y,1)/sqrt(length(y));
cc=corrcoef(xx*b,y);
cc=cc(1,2);

if  plot_results
    figure;
    subplot(1,2,1);
    
    [mask2,rows,cols]=crop_zero(mask);
    imshow(ph(rows,cols),[]);
    
    subplot(1,2,2);
    plot(x,y,'k.');
    hold on;
    plot(x,xx*b,'r-');
    title(sprintf('M = %4.3e',dchi_a2));
     xlabel('C\timescos(2\theta)');
      ylabel('Phase (radian)');
     
end
