if 1
tic;
fname ='meas_MID59_tse_vfl_pssCSspherePol8_FID9136.dat';

%fname ='meas_MID144_tse_vfl_pss_CS_FID8009.dat';
fd224=[];
p=[];
s=[];
i=0;
while 1

%[d,p,s]=reads_fid(fname,-1,'coilScaling.txt');


[d,ptmp,stmp]=readMeasDat(fullfile('f:',fname),'coilScaling.txt',10000,10000*i);

d=(d(1:2:end,:)+d(2:2:end,:))/2;

fd=ifft1c(d,1);

p=[p,ptmp];
s=[s,stmp];


fd224=cat(2,fd224,fd(224,:));

i=i+1;
if size(d,2)<10000
    break;
end

end


p=p(1:32:end);
s=s(1:32:end);
%%
df=single(zeros(max(p)-min(p)+1,max(s)+1,32));
m=zeros(max(p)-min(p)+1,max(s)+1);
for i=1:length(s)
  df(p(i)-min(p)+1,s(i)+1,:)=fd224((i-1)*32+1:i*32);
  m(p(i)-min(p)+1,s(i)+1)=1;
end

end
%%

%[rows,cols, dcomp] = getCalibSizeTSE(m);
cols=121:136;

rows=96:118;

%cols=cols(ceil(end/2)-11:ceil(end/2)+12);


sz=size(df);

wavWeight = 0.002;
toc
tic;

nneg=min(rows)-1+floor(length(rows)/2);
npos=size(m,1)-max(rows)+floor((length(rows)-1)/2);

nskip=npos-nneg+1;
if nskip<=2
 im_pocsspirit=single(zeros(sz(1:3)));
else
 sz2=sz;
sz2(1)=sz(1)+nskip;

  im_pocsspirit=single(zeros(sz2(1:2)));

end
%matlabpool open 64

%%%
%parfor iro=224%1:size(df,1)


DATA =squeeze(df(:,:,:));
pe = size(DATA,2); fe = size(DATA,1); coils = size(DATA,3);

DATAcomp = DATA;%.*repmat(dcomp,[1,1,coils]);
scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;

DATA = DATA/scale_fctr;
DATAcomp = DATAcomp/scale_fctr;

im_dc = ifft2c(DATAcomp);

CalibTyk = 0.01;
kSize = [7,7];


kCalib = DATA(rows,cols,:);
kernel = zeros([kSize,coils,coils]);

[AtA] = corrMatrix(kCalib,kSize);

for n=1:coils
	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
end

GOP = SPIRiT(kernel, 'fft',[fe,pe]);


disp('performing pocs reconstruction');


nIter = 20;

[res_pocs] = pocsSPIRiT(DATA,GOP,nIter,DATA,wavWeight,0);

if nskip<=2
 imtmp= ifft2c(res_pocs);
else
 imtmp=partialFT(permute(res_pocs,[2,1,3]),nneg);
imtmp=permute(imtmp,[2,1,3]);
end

im_pocsspirit=sqrt(sum(abs(imtmp).^2,3));
disp(iro);


%end


%%
matlabpool close
save(sprintf('%s.mat',strtok(fname,'.')),'im_pocsspirit');





