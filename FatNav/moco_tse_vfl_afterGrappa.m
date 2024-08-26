function moco_tse_vfl_afterGrappa(Data,fpro,dfile,nblock,prefix_out)
% prefix: prefix of the rawdata 
% 12/31/2019: scale back to the original voxel size; 12/31/2019
%
% 9/28/2020: replace prefix/prefix with fpro
if ~exist('debug','var')
    debug=false;
end


fov(1)=readsPar(fpro,'asSlice[0].dReadoutFOV');
fov(3)=readsPar(fpro,'asSlice[0].dThickness');
fov(2)=readsPar(fpro,'asSlice[0].dPhaseFOV');


lBaseResolution = readsPar(fpro,'lBaseResolution');
lPhaseEncodingLines =readsPar(fpro,'lPhaseEncodingLines');
lPartitions = readsPar(fpro,'lPartitions');
lImagesPerSlab = readsPar(fpro,'lImagesPerSlab');

lMatrix=[lBaseResolution,lPhaseEncodingLines,lImagesPerSlab,lPartitions];

pnav_sag = readsPar(fpro,'sCuboid.sPosition.dSag');
pnav_cor= readsPar(fpro,'sCuboid.sPosition.dCor');
pnav_tra=readsPar(fpro,'sCuboid.sPosition.dTra');
if isempty(pnav_sag)
    pnav_sag=0;
end
if isempty(pnav_cor)
    pnav_cor=0;
end
if isempty(pnav_tra)
    pnav_tra=0;
end

p_sag = readsPar(fpro,'asSlice[0].sPosition.dSag');
p_cor= readsPar(fpro,'asSlice[0].sPosition.dCor');
p_tra=readsPar(fpro,'asSlice[0].sPosition.dTra');


if isempty(p_sag)
    p_sag=0;
end
if isempty(p_cor)
    p_cor=0;
end
if isempty(p_tra)
    p_tra=0;
end

if isa(Data,'char')
 Data=ri(Data);
end

Data = fft1c(fft1c(Data,3),1);   % PAR to k space; PE already in k space
Nc=32;
fov_os=[fov(1:2),fov(3)*lMatrix(4)/lMatrix(3)];  % fov_os after oversampling along Partition direction
ro_pe_par=[1,-2,-3];
ds = [p_sag,p_cor,p_tra]-[pnav_sag,pnav_cor,pnav_tra];
imSize=lMatrix([1,2,4]);

res=zeros(imSize,'single');

[k2,Data]=get_k_data(Data,dfile,fov_os,ro_pe_par,ds);
        
 nufft=NUFFT3D(k2,1,[0,0,0],imSize,nblock);
 s=tic;
 for i=1:size(Data,2) 
 
  res=res+abs(nufft'*Data(:,i)).^2;
  time_left(i,size(Data,2),toc(s));
 end
 
 
%dfile2=filename(dfile);
%grappa2=filename(grappa);

%prefix_out=[strtok(grappa2,'.'),'_',strtok(dfile2,'.')];
%res=single(abs(res));
res=sqrt(res);
mat2nii_TSE_recon(res,fpro,prefix_out);
%save(prefix_out,'res');


function [k2,Data]=get_k_data(Data,dfile,fov,ro_pe_par,ds)
% fov [3]: field of v for [ro, pe,par]
%
% ro_pe_par: [3]; first,second, third elements for ro, pe, par, respectively,
                 % 1 - x (LR); 2 - y (AP); 3 - z (IS);
% if dfile is empty, the k data will not be transformed.


  
Nc=size(Data,4);
sz=size(Data);
nro=size(Data,1);

Data=reshape(Data,[nro,sz(2)*sz(3),Nc]);


Line0=floor(sz(2)/2);
Partition0=floor(sz(3)/2);
voxsize=fov./sz(1:3);

kpe=((1:sz(2))-Line0)/sz(2)*voxsize(1)/voxsize(2);  % RL
kro=(-nro/2:nro/2-1)/nro;   %AP
kpar=((1:sz(3))-Partition0)/sz(3)*voxsize(1)/voxsize(3); % IS
kpar=reshape(kpar,[1,1,sz(3)]);
kpe=reshape(kpe,[1,sz(2),1]);

kpe=repmat(kpe,[nro,1,sz(3)]);
kpar=repmat(kpar,[nro,sz(2),1]);
kro=repmat(kro',[1,sz(2),sz(3)]);

k=[kro(:),kpe(:),kpar(:)];

for i=1:3
    if ro_pe_par(i)<0
        k(:,i)=-k(:,i);
    end
end

k(:,abs(ro_pe_par))=k;
%imSize(ro_pe_par)=imSize;
voxsize(abs(ro_pe_par))=voxsize;

%%

    if strcmp(dfile(end-3:end),'.mat')% data from afni_motionPar_newBase;
        tmp=load(dfile);
        dfile=[tmp.an',squeeze(tmp.shift([3,1,2],1,:))'];  
        dfile = [(1:size(dfile,1))',dfile]; %to be consistent with .1D file.
    end
    
    
    xform=afni_motionPar2Mat(dfile);
    xform(end+1:sz(3),:)=0; % no motion info for data omitted due to partial Fourier
   
    xform=reshape(xform,[size(xform,1),4,3]);
    xform=permute(xform,[3,2,1]);


    k2=0*k;
    npe=sz(2);
    
    for j=1:size(xform,3)
            ind= (j-1)*npe+1:j*npe;
            ind2= (j-1)*npe*nro+1:j*npe*nro;

       [k2(ind2,:),Data(:,ind,:)] = kspace_xform(k(ind2,:),Data(:,ind,:),xform(:,:,j),voxsize,ds);    
    end    

k2(:,2)=k2(:,2)*voxsize(2)/voxsize(1); %scale back to the original voxel size; 12/31/2019
k2(:,3)=k2(:,3)*voxsize(3)/voxsize(1); %scale back to the original voxel size

%%

Data = reshape(Data,[length(Data(:))/Nc,Nc]);



function [k2,d3] = kspace_xform(k,data,xform,voxsize,ds)
% phase encoding is the 2nd dimension
% k: (nro*npe)*3; -0.5 0.5
% data: nro*npe*nch
% xform: 3*4
% fov: 1*3


%k2=xform(:,1:3)*k';
%k2=k2';

k2=xform(:,1:3)*k';

k2=k2';


k3=reshape(k2',[3,size(data,1),size(data,2)]);
k3=repmat(k3,[1,1,1,size(data,3)]);

k3=permute(k3,[2,3,4,1]);
shift=xform(:,4)'; % phase shift
%shift=shift-((xform(:,1:3)-eye(3))*ds(:))';  % result is worse with - sign
shift=shift+((xform(:,1:3)-eye(3))*ds(:))';

shift=shift./voxsize(1)*2*pi;  %changed from voxsize to voxsize(1); 8/20/2020

shift=reshape(shift,[1,1,1,3]);


shift=repmat(shift,[size(data),1]);


d3 = data.*exp(-1i*sum(shift.*k3,4));  % 


