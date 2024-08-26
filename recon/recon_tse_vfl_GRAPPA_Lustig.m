function recon_tse_vfl(fname,ine)
% fname: *.mat file 
%ine: do grappa on the plane perpendicular to dimensions ine;

prefix=strtok(fname,'.');

lMatrix(1)=readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
lMatrix(2)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
lMatrix(3)=readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');
lMatrix(4)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
Nc=nChan_SiemensProt(fullfile(prefix,[prefix,'.pro']));

load(fname);

if ine==3
 Data=fft1c(Data,1);  % go back to the k-space
end

nro=size(Data,1);

imSize=[nro,lMatrix(2),lMatrix(4)];
[Data3,mask]=sort_kData(Data,Line,Partition,Nc,imSize);
clear Data;

if ine==3
 for i=1:size(Data3,4)
  Data3(:,:,:,i)=ifft1c(Data3(:,:,:,i),3);
  disp(i);
 end
end



p=Line(1:Nc:end)-min(Line);
nskip=max(find(diff(double(p))==0))+1;
p2 = p(nskip+1:end);  % this needs to be changed.

acs_pe=acsLines(unique(p2)+1);

res=zeros(imSize);

for i=1:size(Data3,ine)
    disp(i);
    tic;
    if ine==3
        kCalib=squeeze(Data3(end/2-19:end/2+20,acs_pe,i,:));       
        restmp = GRAPPA(squeeze(Data3(:,:,i,:)),kCalib,[7,7],0);
        res(:,:,par_array(i))=sos(ifft2c(restmp),3);
        
    elseif ine==1
        
        
        kCalib=squeeze(Data3(i,acs_pe,end/2-11:end/2+12,:));    
        dtmp=squeeze(Data3(i,:,:,:));
        dtmp=permute(dtmp,[2,1,3]);
        kCalib=permute(kCalib,[2,1,3]);
        restmp = GRAPPA(dtmp,kCalib,[7,7],0);
    
        restmp=permute(restmp,[2,1,3]);
        res(i,:,:)=sos(ifft2c(restmp),3);
        
    end
    time_left(i,size(Data3,ine),toc);
    
end

fname=unique_name(sprintf('%s_GRAPPA_ine%d.mat',prefix,ine));
save(fname, 'res','-v7.3');

        
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



