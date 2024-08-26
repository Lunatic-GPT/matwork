function [m,DAC,coil]=fm_readcalib

if exist('/home/xiaopeng/labhome/vnmrsys/fastmap','dir')
    
file = '/home/xiaopeng/labhome/vnmrsys/fastmap/fastmap.calib';
else
    file='/data/users/xiaopeng/vnmrsys/fastmap/fastmap.calib';
end

data=textread(file,'%s');
data=reshape(data,[10,7,8]);  %first index measured value, 2-DAC values, 3-coil DAC was applied.
coil=squeeze(data(1,1,:))'; 

DAC=zeros(7,8);
for i=1:7
    for j=1:8
        
    DAC(i,j)=squeeze(str2double(data{2,i,j}));
    end
end

m= zeros(8,7,8);
for i=1:8
    for j=1:7
        for k=1:8
          m(i,j,k)=str2double(data{2+i,j,k});
        end
    end
end



