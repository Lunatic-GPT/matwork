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
