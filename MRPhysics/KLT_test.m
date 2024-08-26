nts=1000;
N=100;
c=0.7;
ref=load('refg_10_20_80_10cy_TR1.0');
ref=ref(1:1000);
ref=ref/max(ref);
ts=zeros(nts,N);
for j=1:nts
ts(j,1)=randn_white(1)/(1-c^2);
 for i=1:N-1
    ts(j,i+1)=c*ts(j,i)+randn_white(1);
 end
end


fts=fft(ts,[],2);

figure;plot(abs(fts(5,:)));

disp(var(ts)*(1-c^2));
disp(var(randn_white(N)));  %they should be the same

sig=zeros(N,N);

for i=1:N
    for j=1:N
        sig(i,j)=1/(1-c^2)*c^(abs(i-j));
    end
end

[V,mat]=KLT_matrix(ts);

[U,D]=eig(sig);

ts2=U'*(ts');
  
figure;subplot(1,2,1);imshow(ts',[]);
subplot(1,2,2);imshow(abs(ts2),[]);
xlim([0,110]);