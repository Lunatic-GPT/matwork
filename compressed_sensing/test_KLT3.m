
in=sum(data,4);
nm=sum(mask,4);

m2=(nm==size(mask,4));
d3=zeros(length(find(m2>0)),152);
for i=1:152
    tmp=data(:,:,:,i);
    d3(:,i)=tmp(m2);
end

d3orig=d3;
%d3=abs(d3);
[v,mat]=KLT_matrix(d3);
d3m=mean(d3,2);
d3m=repmat(d3m,[1,152]);
b=v'*transpose(d3-d3m);
b2=fft(d3-d3m,[],2);
b2=b2';
figure;imshow(abs(b(:,1:300)),[]);

figure;imshow(abs(b2(2:end,1:300)),[]);

thr=max(abs(b2(:)))*0.05;
frac=length(find(abs(b2(:))>thr))/length(b2(:))

thr=max(abs(b(:)))*0.05;
frac=length(find(abs(b(:))>thr))/length(b2(:))




