function kdata=recon_tse_vfl_Siemens(fname,prefix_out,fpro,debug)
% fname: *.mat file
% use the Siemens grappa setting, i.e., do IFT on both fully sampled
% dimensions before GRAPPA.
% 9/28/2020: add fpro parameter, was using prefix/prefix.pro
if ~exist('debug','var')
    debug=false;
end



lMatrix(1)=readsPar(fpro,'lBaseResolution');
lMatrix(2)=readsPar(fpro,'lPhaseEncodingLines');
lMatrix(3)=readsPar(fpro,'lImagesPerSlab');
lMatrix(4)=readsPar(fpro,'lPartitions');
Nc=nChan_SiemensProt(fpro);

load(fname);

if debug
    Data=Data(end/2:end/2+1,:,:,:);
end


nro=size(Data,1);

imSize=[nro,lMatrix(2),lMatrix(4)];

[Data3,Slines]=sort_kData(Data,Line,Partition,Nc,imSize);    

clear Data;

for i=1:size(Data3,4)
   Data3(:,:,:,i)=ifft1c(Data3(:,:,:,i),3);
end

[res,kdata,coef]=recon_grappa2D_Siemens(Data3,Slines,imSize(2));

%coef_ch0=reshape(squeeze(coef(:,:,:,1,:)),[size(coef,1),Nc,size(coef,2)/Nc,size(coef,3),size(coef,5)]);
%save(fname, 'res','coef','-v7.3');% coef takes a lot of memory
save([prefix_out,'.mat'], 'res','-v7.3');
% 
% for i=1:size(kdata,4)
%    fname=sprintf('%s_GRAPPA_Siem_Ch%d.mat',prefix,i); 
%    kdata_sc=kdata(:,:,:,i);
%    save(fname,'kdata_sc');  
% end


function [Data3,Slines,dnoise]=sort_kData(Data2,Line,Partition,Nc,imSize)


%p=unique(Line);

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);
n=size(Data2,2);
n=ceil(n/32)*32;

Data2(:,end+1:n)=0;
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



