function recon_ESPIRiT_3DData(fname,iro_array,nmaps,lMatrix,eigThresh_k)
%fname='meas_MID85_tse_vfl_pss_FatNav_FID23810.mat';

tcg=tic;

if ~exist('lMatrix','var')
lMatrix = [512,405,244,256];
end

if ~exist('eigThresh_k','var')
    eigThresh_k=0.02;
end

load(fname);
nro=size(Data,1);
Nc=32;

prefix=strtok(fname,'.');

imSize=[lMatrix(2),nro,lMatrix(4)];

% only for this time; below
res=zeros([lMatrix(2),length(iro_array),lMatrix(4),nmaps]);
scl=zeros(1,length(iro_array));
for i=1:length(iro_array)
iro=iro_array(i);
[Data3,mask,DCalib]=sort_kData(Data(iro,:),Line,Partition,Nc,imSize);

mask=repmat(mask,[1,1,Nc]);
FT=p2DFT(mask,size(mask),1,2);

scl(i)=max(abs(Data3(:)));
Data3=Data3/scl(i)+eps;
DCalib=DCalib/scl(i);  %this scaling does not matter
kSize=[7,7];
eigThresh_im = 0.99;
%eigThresh_k = 0.02;  %0.05 is not good, use 0.02;

    [kernel,s] = dat2Kernel(DCalib,kSize);
    idx = max(find(s >= s(1)*eigThresh_k));
    if idx>200% most likely no signal; just noise
       % continue;
    end

    [M,W] = kernelEig(kernel(:,:,:,1:idx),imSize([1,3]));
    maps = M(:,:,:,end-nmaps+1:end);
    weights = double((W(:,:,end-nmaps+1:end) >  eigThresh_im));
   % weights = ones(size(W(:,:,end-nmaps+1:end)));
    save Data3_tmp maps weights Data3 mask 
    disp('Construct ESP ...');   
    tic;
    ESP = ESPIRiT(maps,weights);  
    disp('ESP construction done');
    toc;

XOP = Wavelet('Daubechies',4,6);
disp('Start CGL1ESPIRiT');
tcg=tic;
nIterCG = 5;
im0=ESP'*(FT'*Data3);
 %  res(:,i,:,:)= cgL1ESPIRiT(double(Data3), zeros(imSize(1),imSize(3),nmaps), FT,ESP, nIterCG,XOP,0,0,1);

   res(:,i,:,:)= cgL1ESPIRiT(double(Data3), im0, FT,ESP, nIterCG,XOP,0,0,1);


   res(:,i,:,:)=res(:,i,:,:)*scl(i);
end
%tmp = cgL1ESPIRiT(double(DATA), zeros(sx,sy,2), FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);


disp('CGL1ESPIRiT finished');
toc(tcg);


fname=unique_name(sprintf('%s_ESPIRiT_%dmaps_%d_%d_%s.mat',prefix,nmaps,iro_array(1),iro_array(end),num2str(eigThresh_k)));

save(fname, 'res','maps','W','-v7.3');


function [Data3,mask,DCalib]=sort_kData(Data2,Line,Partition,Nc,imSize)

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);



Data2=reshape(Data2,[Nc,length(Data2)/Nc]);

nskip=max(find(diff(double(p))==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

Data3 = zeros([imSize(1),imSize(3),Nc],'single');
mask=zeros([imSize(1),imSize(3)]);
%%

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

for i=1:length(s2)
    Data3(p2(i)+1,s2(i)+1,:)=Data2(:,nskip+i);
    mask(p2(i)+1,s2(i)+1,:) = 1;
end


acs_pe=acsLines(unique(p2)+1);
acs_par = imSize(3)/2-25:imSize(3)/2+26;


DCalib=Data3(acs_pe,acs_par,:);




