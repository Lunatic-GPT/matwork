function c=centerFinderPhase(cd,range,i_rad_roi,i_rad_bg)
% cd: complex data
% rad: radius of circle mask or width of a square (2*rad+1).
% range: the size of the square around center to search the center is
% (2*range+1)*(2*range+1)
% if circ_square = 0: use circle otherwise use a square.
% this method does not work well

cd=permute(cd,[2,1,3]);




c0=ceil((size(cd)+1)/2);

rl=zeros(2*range+1,2*range+1);

for i=-range:range
    for j=-range:range
        
        c=[c0(1)+i,c0(2)+j];
        m_vessel=mask_circle(size(cd),i_rad_roi,c(1:2),1);
        
        m_large = mask_circle(size(cd),i_rad_bg,c(1:2),1);
        
      %  m_large=permute(m_large,[2,1,3]);
      %  m_vessel=permute(m_vessel,[2,1,3]);
        
        rl(i+range+1,j+range+1)= magMoment_outSidePattern(cd,m_large>0&m_vessel==0,c0+[i,j],false);
        
    end
end


[tmp,ind_min]=min(real(rl(:)));

cshift=ind2subb(size(rl),ind_min);

c=c0+cshift+[-range-1,-range-1];
m_vessel=mask_circle(size(cd),i_rad_roi,c(1:2),1);
m_large = mask_circle(size(cd),i_rad_bg,c(1:2),1);
res=magMoment_outSidePattern(cd,m_large>0&m_vessel==0,c,true);



function resid=magMoment_outSidePattern(ph,mask,center,plot_results)

% B0 projection is along the second dimension
% center: indices for the center
% dchi_a2 unit: ppm*mm2
% ph in degree
%
ph=permute(ph,[2,1,3]);
mask=permute(mask,[2,1,3]);

pos=find(mask>0);
pos=ind2subb(size(mask),pos);
ndim=ndims(mask);
pos=pos-repmat(center,size(pos,1),1);
%pos=pos.*repmat(voxSize(1:ndim),size(pos,1),1);

phi=atan2(pos(:,1),pos(:,2));

x=cos(2*phi)./sum(pos(:,1:2).^2,2);

y=ph(mask>0);



xx=[x(:),ones(length(x),1)];

b=xx\y;
resid=sos(xx*b-y,1)/sqrt(length(y));


if plot_results
    figure;plot(x,y,'k.');
    hold on;
    plot(x,xx*b,'r-');
    
end
