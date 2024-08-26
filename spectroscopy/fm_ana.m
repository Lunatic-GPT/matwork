a=read_fid('fid');
gfs=readPar('procpar','gfs');
gf=readPar('procpar','gf');
sw=readPar('procpar','sw');
vox1=readPar('procpar','vox1');%y
vox2=readPar('procpar','vox2');%x
vox3=readPar('procpar','vox3');%z
vox1=vox1/10;
vox2=vox2/10;
vox3=vox3/10;
pos(1)=readPar('procpar','pos1');%y
pos(2)=readPar('procpar','pos2');%x
pos(3)=readPar('procpar','pos3');%z
tau=readPar('procpar','tau');%z

lro=readPar('procpar','lro');
a=squeeze(a);

i=0:size(a,1)-1;
g=exp(-(i/sw-gfs).^2/gf^2);

ag=a.*repmat(g',[1,size(a,2)]);
ags=fftshift(ag,1);
fa=fft(ags,[],1);
fa=fftshift(fa,1);

figure;
for i=1:6
 subplot(2,3,i);
 plot(real(a(:,2*i-1)),'r-');
 hold on;plot(real(a(:,2*i)),'b-');
 yl=ylim;
 hold on;plot(g*yl(2)/2,'k--');
 
end
  
figure;
for i=1:6
 subplot(2,3,i);
 plot(real(fa(:,2*i-1)),'r-');
 hold on;plot(real(fa(:,2*i)),'b-');
end

figure;
for i=1:6
 subplot(2,3,i);
 plot(phase(fa(:,2*i-1)),'r-');
 hold on;plot(phase(fa(:,2*i)),'b-');
end


%% calculate the shim fields.

change=read_shimchanges('tao');
[m,cc]=fm_readcorr;
m(abs(cc)<0.93)=0;
k=m*change';
disp(k);
mk2tok1=matrix_offset_k2tok1([-pos(2),-pos(1),pos(3)]);
ka=mk2tok1*k(4:end);
k(1:3)=ka+k(1:3);


theta=[pi/2,pi/2,pi/4,pi/4,pi/4,pi/4];
phi=[pi/4,-pi/4,pi/2,-pi/2,0,pi];

dr=linspace(-1,1,100);
Bshim = zeros(100,6);
for i=1:6
    for j=1:length(dr)
       Bshim(j,i)=fm_shimfield(theta(i),phi(i),k,pos,dr(j))*tau*2*pi;
    end
end

r(1)=sqrt(vox1^2+vox2^2);
r(2)=sqrt(vox1^2+vox2^2);
r(3)=sqrt(vox1^2+vox3^2);
r(4)=sqrt(vox1^2+vox3^2);
r(5)=sqrt(vox2^2+vox3^2);
r(6)=sqrt(vox2^2+vox3^2);

ttl={'X Y','X -Y','Y Z','-Y Z','X Z','-X Z'};

r_array=linspace(-lro/2,lro/2,size(fa,1));
figure;
for i=1:6
 subplot(2,3,i);
 y = phase(fa(:,2*i-1))-phase(fa(:,2*i));
 plot(r_array,y,'r-');
 hold on;plot([-r(i)/2,-r(i)/2],ylim,'k--');
 plot([r(i)/2,r(i)/2],ylim,'k--');
 title(ttl{i});
 plot(fliplr(dr),Bshim(:,i)+y(size(fa,1)/2),'k');
end

%%








