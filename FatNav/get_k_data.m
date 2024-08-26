function [k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov,lMatrix,ro_pe_par,ds)
% fov [3]: field of v for [ro, pe,par]
% lMatrix: matrix size for [pe,par];
%
% ro_pe_par: [3]; first,second, third elements for ro, pe, par, respectively,
                 % 1 - x (LR); 2 - y (AP); 3 - z (IS);
% if dfile is empty, the k data will not be transformed.


if ~exist('ro_pe_par','var')
    ro_pe_par=[2,1,3];
end

Line=double(Line);
Partition=double(Partition);


nkeep=floor(length(Line)/Nc)*Nc;
Line=Line(1:nkeep);
Partition=Partition(1:nkeep);
Data=Data(:,1:nkeep);

%xformi=invert_m(xform);  %from image to base 3*4*n

  
p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);
nro=size(Data,1);




nskip=max(find(diff(p)==0))+1;% last index for phase correction scan
if isempty(nskip)% no phase correction scan; only skip the GRAPPA noise scan
    if abs(p(2)-p(1))>1     
        if length(unique(s))*length(unique(p))*Nc ==size(Data,2) % no noise scan
            nskip=0;
        else      
            nskip=1;
        end
    else       
        nskip=0;
    end
end

Data=reshape(Data,[nro,Nc,size(Data,2)/Nc]);
Data=permute(Data,[1,3,2]);
%%

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);
Data = Data(:,nskip+1:end,:);

m=zeros(max(p2)+1,max(s2)+1);


for i=1:length(s2)

  m(p2(i)+1,s2(i)+1)=1;
end

%         
%   d2=zeros(nro,max(p2)+1,max(s2)+1,'single');    
%         for i=1:length(s2)
%             d2(:,p2(i)+1,s2(i)+1)=sos(Data(:,nskip+i,:),3);
%         end
%  %
% ips0=max_ind(squeeze(sum(d2,1)));
% Line0=ips0(1)-1; % PE = 0 ;
% Partition0=ips0(2)-1;  % par = 0;
Line0=floor(lMatrix(1)/2);
Partition0=floor(lMatrix(2)/2);


    
%ro_center=nro/2-25:nro/2+26;
voxsize=fov./[nro,lMatrix];
kpe=(p2-Line0)/lMatrix(1)*voxsize(1)/voxsize(2);  % RL
kro=(-nro/2:nro/2-1)/nro;   %AP
kpar=(s2-Partition0)/lMatrix(2)*voxsize(1)/voxsize(3); % IS

%imSize=[nro,lMatrix];

%voxsize=fov([2,1,3])./[lMatrix(1),nro,lMatrix(2)];
%imSize=[lMatrix(1),nro,lMatrix(2)];

kpe=reshape(kpe,[1,length(kpe)]);
kpar=reshape(kpar,[1,length(kpar)]);
kpe=repmat(kpe,[nro,1]);
kpar=repmat(kpar,[nro,1]);
kro=repmat(kro',[1,size(kpe,2)]);

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

if ~isempty(dfile)
    if strcmp(dfile(end-3:end),'.mat')% data from afni_motionPar_newBase;
        tmp=load(dfile);
        dfile=[tmp.an',squeeze(tmp.shift([3,1,2],1,:))'];  
        dfile = [(1:size(dfile,1))',dfile]; %to be consistent with .1D file.
    end
    
    xform=afni_motionPar2Mat(dfile);  % should get the same value with the following command
    xform=reshape(xform,[size(xform,1),4,3]);
    xform=permute(xform,[3,2,1]);
else
    xform=repmat([1,0,0,0;0,1,0,0;0,0,1,0],[1,1,length(unique(s))]);
end

if size(xform,3)~=length(unique(s))
    error('size(xform,3)~=length(unique(s))');
end

npe = length(unique(p2));

k2=0*k;

    for j=1:size(xform,3)
        if j<size(xform,3)
            ind= (j-1)*npe+1:j*npe;
            ind2= (j-1)*npe*nro+1:j*npe*nro;
        else
            ind= (j-1)*npe+1:size(Data,2);  % some data are lost for the last TR
            ind2= (j-1)*npe*nro+1:size(k2,1);
        end
       [k2(ind2,:),Data(:,ind,:)] = kspace_xform(k(ind2,:),Data(:,ind,:),xform(:,:,j),voxsize,ds);    
    end    


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

shift=shift./voxsize*2*pi;

shift=reshape(shift,[1,1,1,3]);


shift=repmat(shift,[size(data),1]);


d3 = data.*exp(-1i*sum(shift.*k3,4));  % 

