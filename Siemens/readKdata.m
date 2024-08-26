fname='meas_MID48_tse_vfl_pssCSsphere_FID8552.dat';

pall=[];
sall=[];
i=0;
n=32*1000;
ix=43:415;
fdall=[];
while 1
[d,p,s]=readMeasDat(fname,'coilScaling.txt',n,i*n);

d=d(1:2:end,:)+d(2:2:end,:);

fd=ifft1c(d,1);
fdall=cat(2,fdall,fd(ix,:));

pall=[pall;p(:)];
sall=[sall;s(:)];

i=i+1;
if isempty(d)|| size(d,2)<n
    break;
end

end
%save fd224_8159 fd224 pall sall;

fdall=single(fdall);
pall=pall(1:32:end);
sall=sall(1:32:end);

save(strtok(fname,'.'),'pall','sall','fdall','-v7.3');
