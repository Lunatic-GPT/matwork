function recon_tse_vfl(fname,ine,debug)
% fname: *.mat file
%ine: do grappa on the plane perpendicular to dimensions ine;
if ~exist('debug','var')
    debug=false;
end


prefix=strtok(fname,'.');

lMatrix(1)=readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
lMatrix(2)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
lMatrix(3)=readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');
lMatrix(4)=readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
Nc=nChan_SiemensProt(fullfile(prefix,[prefix,'.pro']));

load(fname);

if debug
  ine=1;
  Data=Data(end/2,:,:,:);
  Dnoise=Data(:,1:Nc);
end

if ine==3
    Data=fft1c(Data,1);  % go back to the k-space
end

nro=size(Data,1);

imSize=[nro,lMatrix(2),lMatrix(4)];
[Data3,Slines]=sort_kData(Data,Line,Partition,Nc,imSize);
clear Data;

if ine==3
    for i=1:size(Data3,4)
        Data3(:,:,:,i)=ifft1c(Data3(:,:,:,i),3);
    end
end

if ine==3
    [res,kdata]=recon_grappa2D(Data3,Slines,imSize(2));
elseif ine==1
    [res,kdata]=recon_grappa2D(permute(Data3,[3,2,1,4]),Slines,imSize(2));
    res=permute(res,[3,2,1,4]);
end


fname=unique_name(sprintf('%s_GRAPPA_ine%d.mat',prefix,ine));
save(fname, 'res','-v7.3');

for i=1:size(kdata,4)
   fname=sprintf('%s_GRAPPA_Ch%d.mat',prefix,i); 
   kdata_sc=kdata(:,:,:,i);
   save(fname,'kdata_sc');  
end


function [Data3,Slines,dnoise]=sort_kData(Data2,Line,Partition,Nc,imSize)


%p=unique(Line);

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);

Data2=reshape(Data2,[size(Data2,1),Nc,size(Data2,2)/Nc]);

nskip=max(find(diff(double(p))==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

dnoise=Data2(:,:,1);

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

p2_unq=unique(p2);

Data3 = zeros([imSize(1),length(p2_unq),imSize(3),Nc],'single');

for i=1:length(s2)
    Data3(:,p2(i)==p2_unq,s2(i)+1,:)=Data2(:,:,nskip+i);
end

Slines=p2_unq+min(Line)+1;


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



