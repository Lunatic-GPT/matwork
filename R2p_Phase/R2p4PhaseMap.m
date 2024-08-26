function R2p=R2p4PhaseMap(ph,TE,usFactor,mask)

% ph: in units of radian; 
% usFactor: undersampling factors in x, y, and z directions

R2p=0*ph;
gamma=4258*2*pi;

sz=size(ph);
for i=1:size(ph,1)
    disp(i);
    for j=1:size(ph,2)
        for k=1:size(ph,3)
            
            i1=[i,j,k]-floor(usFactor/2);
            i2=[i,j,k]+floor(usFactor/2);
            
            i1(i1<1)=1;
            i2(i2>sz)=sz(i2>sz);
            
            sel=mask(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3));
            
            if sum(sel(:))>=5
                
                tmp=ph(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3));
                sd=std(tmp(sel>0))/TE/gamma;
                R2p(i,j,k)=BzSD2R2p(sd);
            end

        end
    end
end



