function recon_ESPIRiT_3DData_PaFT(fname,par_array,nmaps,lMatrix,eigThresh_k)
% perform ESPIRiT along ro*pe, zero filled IFT along par encoding

prefix=strtok(fname,'.');
lMatrix(1)=readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
lMatrix(2)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
lMatrix(3)=readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');
lMatrix(4)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
Nc=nChan_SiemensProt(fullfile(prefix,[prefix,'.pro']));




tcg=tic;


if ~exist('eigThresh_k','var')
    eigThresh_k=0.01;
end

load(fname);
nro=size(Data,1);

imSize=[nro,lMatrix(2),lMatrix(4)];


% only for this time; below

scl=zeros(1,length(par_array));
[Data3,mask]=sort_kData(Data,Line,Partition,Nc,imSize);
for i=1:size(Data3,4)
Data3(:,:,:,i)=ifft1c(Data3(:,:,:,i),3);
disp(i);
end

p=Line(1:Nc:end)-min(Line);
nskip=max(find(diff(double(p))==0))+1;
p2 = p(nskip+1:end);  % this needs to be changed.

    mask=repmat(mask(:,1)',[nro,1,Nc]);
    
for i=1:length(par_array)
%%
    
    Data4=squeeze(Data3(:,:,par_array(i),:));
    
  %  Data4=do_hanning(Data4,2);
   
    Data4=fft1c(Data4,1);
    
    scl(i)=max(abs(Data4(:)));

    Data4=Data4/scl(i)+eps;
    
    

kSize=[5,5];
eigThresh_im = 0.98;
%eigThresh_k = 0.02;  %0.05 is not good, use 0.02;
acs_pe=acsLines(unique(p2)+1);

DCalib=Data4(end/2-19:end/2+20,acs_pe,:);

    [kernel,s] = dat2Kernel(DCalib,kSize);
    
     eigThresh_k=0.01;
    idx = max(find(s >= s(1)*eigThresh_k));
    if idx>200% most likely no signal; just noise
       % continue;
    end
   idx=125;
    imSize2=[512,512,32];
    mask2=crop(mask,imSize2);
    FT=p2DFT(mask2,size(mask2),1,2);
    Data5=crop(Data4,imSize2);
      
    [M,W] = kernelEig(kernel(:,:,:,1:idx),imSize2([1,2]));
 
    maps = M(:,:,:,end-nmaps+1:end);
    weights = double((W(:,:,end-nmaps+1:end) >  eigThresh_im));
   % weights = ones(size(W(:,:,end-nmaps+1:end)));
 %   save Data3_tmp maps weights Data3 mask 
    disp('Construct ESP ...');   
    tic;
    ESP = ESPIRiT(maps,weights);  
    disp('ESP construction done');
    toc;

XOP = Wavelet('Daubechies',4,6);
disp('Start CGL1ESPIRiT');
tcg=tic;
nIterCG = 5;

im0=ESP'*(FT'*squeeze(Data6));
%  res(:,i,:,:)= cgL1ESPIRiT(double(Data3), zeros(imSize(1),imSize(3),nmaps), FT,ESP, nIterCG,XOP,0,0,1);
   
res= cgL1ESPIRiT(double(Data5), im0, FT,ESP, nIterCG,XOP,0,0,3)*scl(i);
im2=ESP*res;
end
%tmp = cgL1ESPIRiT(double(DATA), zeros(sx,sy,2), FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);


disp('CGL1ESPIRiT finished');
toc(tcg);


fname=unique_name(sprintf('%s_ESPIRiT_%dmaps_ksz%d_PaFT_%d_%d_%s.mat',prefix,nmaps,kSize(1),par_array(1),par_array(end),num2str(eigThresh_k)));

save(fname, 'res','maps','W','-v7.3');


function [Data3,mask]=sort_kData(Data2,Line,Partition,Nc,imSize)

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);

Data2=reshape(Data2,[size(Data2,1),Nc,size(Data2,2)/Nc]);

nskip=max(find(diff(double(p))==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

Data3 = zeros([imSize(1),imSize(2),imSize(3),Nc],'single');
mask=zeros([imSize(2),imSize(3)]);
%%

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

for i=1:length(s2)
    Data3(:,p2(i)+1,s2(i)+1,:)=Data2(:,:,nskip+i);
    mask(p2(i)+1,s2(i)+1) = 1;
end



function mask=get_mask(Line,Partition,Nc,imSize)

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);


nskip=max(find(diff(double(p))==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

mask=zeros([imSize(2),imSize(3)]);
%%

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

for i=1:length(s2)
    mask(p2(i)+1,s2(i)+1) = 1;
end



