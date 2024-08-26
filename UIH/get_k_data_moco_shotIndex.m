function [k2,Data]=get_k_data_moco_shotIndex(Data,dfile,shotIndex,fov,ro_pe_par,ds)
% Data: nro*npe*npar*Nch
% fov [3]: field of v for [ro, pe,par]
%
% ro_pe_par: [3]; first,second, third elements for ro, pe, par, respectively,
% 1 - x (RL); 2 - y (AP); 3 - z (IS);
% -1 - -x; -2 - -y; -3 - -z
% if dfile is empty, the k data will not be transformed.

% ds: FOV center difference between main and navigator images
% ds = p_main-p_nav; (dicom convention)

Nc=size(Data,4);
sz=size(Data);
nro=size(Data,1);

Data=reshape(Data,[nro,sz(2)*sz(3),Nc]);

Line0=floor(sz(2)/2);
Partition0=floor(sz(3)/2);
voxsize=fov./sz(1:3);

kpe=((0:sz(2)-1)-Line0)/sz(2)*voxsize(1)/voxsize(2);  % RL
kro=(-nro/2:nro/2-1)/nro;   %AP
kpar=((0:sz(3)-1)-Partition0)/sz(3)*voxsize(1)/voxsize(3); % IS
kpar=reshape(kpar,[1,1,sz(3)]);
kpe=reshape(kpe,[1,sz(2),1]);

kpe=repmat(kpe,[nro,1,sz(3)]);
kpar=repmat(kpar,[nro,sz(2),1]);
kro=repmat(kro',[1,sz(2),sz(3)]);

k=reshape([kro(:),kpe(:),kpar(:)],[nro,sz(2)*sz(3),3]);

for i=1:3
    if ro_pe_par(i)<0
        k(:,:,i)=-k(:,:,i);
    end
end

k(:,:,abs(ro_pe_par))=k;  %[x,y,z]

voxsize(abs(ro_pe_par))=voxsize; %x,y,z

%%
if isempty(dfile)
    dfile = zeros(max(shotIndex(:)),7);
end

xform=afni_motionPar2Mat(dfile);
xform=reshape(xform,[size(xform,1),4,3]);
xform=permute(xform,[3,2,1]);

if size(xform,3)~=max(shotIndex(:))
    error('motion parameters size error')
end

k2=0*k;

shotIndex2=setShotIndex_NN(shotIndex);

for j=1:size(xform,3)
    ind= find(shotIndex2==j);
    [k2(:,ind,:),Data(:,ind,:)] = kspace_xform(k(:,ind,:),Data(:,ind,:),xform(:,:,j),voxsize,ds);
end

k2(:,:,2)=k2(:,:,2)*voxsize(2)/voxsize(1); %scale back to the original voxel size; 12/31/2019
k2(:,:,3)=k2(:,:,3)*voxsize(3)/voxsize(1); %scale back to the original voxel size

%%

Data = reshape(Data,[length(Data(:))/Nc,Nc]);
k2 = reshape(k2,[length(k2(:))/3,3]);



function shotIndex2=setShotIndex_NN(shotIndex)
  
   
sz=size(shotIndex);

shotIndex2=shotIndex;

shotIndex_array=shotIndex(shotIndex>0);
ind=ind2subb(sz,find(shotIndex>0));

for i=1:sz(1)
    for j=1:sz(2)
        
        dsquare=(i-ind(:,1)).^2+(j-ind(:,2)).^2;
        
        [~,ind_min]=min(dsquare);
        
        shotIndex2(i,j)=shotIndex_array(ind_min);

    end
end



function [k2,d3] = kspace_xform(k,data,xform,voxsize,ds)
% phase encoding is the 2nd dimension
% k: nro*npe*3; -0.5 0.5
% data: nro*npe*nch
% xform: 3*4
% fov: 1*3


%k2=xform(:,1:3)*k';
%k2=k2';

sz=size(k);

k=reshape(k,[prod(sz(1:2)),3]);

k2=xform(:,1:3)*k';


k3=reshape(k2,[3,size(data,1),size(data,2)]);
k3=repmat(k3,[1,1,1,size(data,3)]);

k3=permute(k3,[2,3,4,1]);
shift=xform(:,4)'; % phase shift
%shift=shift-((xform(:,1:3)-eye(3))*ds(:))';  % result is worse with - sign
shift=shift+((xform(:,1:3)-eye(3))*ds(:))';

shift=shift./voxsize(1)*2*pi;  %changed from voxsize to voxsize(1); 8/20/2020

shift=reshape(shift,[1,1,1,3]);


shift=repmat(shift,[size(data),1]);


d3 = data.*exp(-1i*sum(shift.*k3,4));  %

k2=reshape(k2',[sz(1:2),3]);

