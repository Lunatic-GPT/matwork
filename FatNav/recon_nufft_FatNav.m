function recon_nufft_FatNav(fname,dfile,fov,lMatrix)
% fname='meas_MID85_tse_vfl_pss_FatNav_FID23810.mat';
% dfile = 'Motion_MID85.1D';
% fov = [210.0,165.7,97.6];
% lMatrix = [512,405,244,256];

% fname='meas_MID18_tse_vfl_pss_FatNav_motion_FID27680.mat';
% dfile='Motion_MID18.1D';
% fov=[256.0,202.0,124.0];
% lMatrix=[256,207,124,128];


load(fname);
Data = fft1c(Data,1);   % to k space;
nro=size(Data,1);
voxsize=fov./[nro,lMatrix(2),lMatrix(3)];

prefix=strtok(fname,'.');
imSize=[lMatrix(2),nro,lMatrix(4)];
 Nc=32;

%%
[k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov,lMatrix);  %
%%
nufft=NUFFT3D(k2,1,[0,0,0],imSize);
    
a=nufft'*Data;

save([prefix,'_nufft.mat'],'a');


function [k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,voxsize,lMatrix)

%%
  dfile =load(dfile);    
  xform=afni_motionPar2Mat(dfile(:,2:7));  % should get the same value with the following command  
  
xform=reshape(xform,[size(xform,1),4,3]);
xform=permute(xform,[3,2,1]);

Line=double(Line);
Partition=double(Partition);
xformi=invert_m(xform);  %from image to base 3*4*n

  
p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);
nro=size(Data,1);



nskip=max(find(diff(p)==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

Data=reshape(Data,[nro,Nc,size(Data,2)/Nc]);
Data=permute(Data,[1,3,2]);


p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

m=zeros(max(p2)+1,max(s2)+1);


for i=1:length(s2)

  m(p2(i)+1,s2(i)+1)=1;
end

        
  d2=zeros(nro,max(p2)+1,max(s2)+1,'single');
        
        
        for i=1:length(s2)
            d2(:,p2(i)+1,s2(i)+1)=sos(Data(:,nskip+i,:),3);
        end
 %
ips0=max_ind(squeeze(sum(d2,1)));
Line0=ips0(1)-1; % PE = 0 ;
Partition0=ips0(2)-1;  % par = 0;


    
%ro_center=nro/2-25:nro/2+26;

kx=(p2-Line0)/lMatrix(2);

ky=(-nro/2:nro/2-1)/nro;

kz=(s2-Partition0)/lMatrix(4);

%



kx=repmat(kx,[nro,1]);
kz=repmat(kz,[nro,1]);
ky=repmat(ky',[1,length(kx)]);


k=[kx(:),ky(:),kz(:)];

voxsize=voxsize([2,1,3]);

if size(xform,3)~=length(unique(s))
    error('number of transformed wrong');
end

ns=length(unique(s2));
npe = length(unique(p2));

k2=0*k;

    for j=1:size(xform,3)
       ind= (j-1)*npe+1:j*npe;
       ind=ind+nskip;
       ind2= (j-1)*npe*nro+1:j*npe*nro;    
       [k2(ind2,:),Data(:,ind,:)] = kspace_xform(k(ind2,:),Data(:,ind,:),xformi(:,:,j),voxsize);    
    end    




Data=Data(:,nskip+1:end,:);
Data = reshape(Data,[length(Data(:))/Nc,Nc]);
