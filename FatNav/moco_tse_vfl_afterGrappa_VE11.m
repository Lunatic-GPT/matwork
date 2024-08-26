function moco_tse_vfl_afterGrappa_VE11(Data,fname,dfile,nblock,prefix_out)
% prefix: prefix of the rawdata
% 12/31/2019: scale back to the original voxel size; 12/31/2019
%
% 9/28/2020: replace prefix/prefix with fpro
% 12/1/2020:
% iRep: same size as the squeeze(Data(1,:,:,1)), the repetitions that the data
% should be assigned to.  range: 1 - number of Nav images/number of rows in
% dfile
% Data is output from recon_tse_vfl_Siemens_VE11


a=mapVBVD(fname);

fov(1)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
fov(2)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
fov(3)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;


lMatrix(1)=a{2}.hdr.Config.NImageCols;
lMatrix(2)=a{2}.hdr.Config.NImageLins;
lMatrix(3)=a{2}.hdr.Config.NoImagesPerSlab;
lMatrix(4)=a{2}.hdr.Config.NoOfFourierPartitions;


ro_pe_par=get_ro_pe_par(a);

Data = fft1c(fft1c(Data,3),1);   % PAR to k space; PE already in k space

fov_os=[fov(1:2),fov(3)*lMatrix(4)/lMatrix(3)];  % fov_os after oversampling along Partition direction

ds = [0,0,0]; %Fatnav and water images have the same FOV
imSize=lMatrix([1,2,4]);

res=zeros(imSize,'single');
iRep=get_iRep(a,lMatrix(2),lMatrix(4));

[k2,Data]=get_k_data(Data,iRep,dfile,ro_pe_par,fov_os,ds);

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
%save(prefix_out,'res');

o=get_rotmat(a,size(res));
mat2niigz(res,'',prefix_out,false,o);


function o=get_rotmat(a,sz)
     fov(1)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
     fov(2)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
     fov(3)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
     
    pos=get_pos(a);
     o.voxsize=fov./sz(:)';
     o.center=pos;
     
     o.rotmat=get_ro_pe_par(a);     
     o.orient='SPR';  
     o.sz=sz;
     o=center2pos_o(o,o.center);

function iRep=get_iRep(a,npe,npar)

lin_par_ref=[a{2}.refscan.Lin(:),a{2}.refscan.Par(:)];
lin_par=[a{2}.image.Lin(:),a{2}.image.Par(:)];
[~,id]=setdiff(lin_par_ref,lin_par,'rows');
ts=cat(2,a{2}.refscan.timestamp(id),a{2}.image.timestamp);
Lin=cat(2,a{2}.refscan.Lin(id),a{2}.image.Lin);
Par=cat(2,a{2}.refscan.Par(id),a{2}.image.Par);

%%
ts_sort=sort(ts);
ind=find(diff(ts_sort)>500)+1;  %the starting index of the next repetition
ind=[1,ind];

ts_start=ts_sort(ind);

%%

iRep=zeros(npe,npar);
tsmap=zeros(npe,npar);

for i=1:npe
    for j=1:npar
        
        dsquare=(i-Lin).^2+(j-Par).^2;
        
        [~,ind]=min(dsquare);
        
        iRep(i,j)=find(ts(ind)-ts_start>=0,1,'last');
        if (min(dsquare)==0)
            tsmap(i,j)=iRep(i,j);
        end
    end
end
%%

function [k2,Data]=get_k_data(Data,iRep,dfile,ro_pe_par,fov,ds)
% fov [3]: field of v for [ro, pe,par]
%
% ro_pe_par: directions for first, second, and third dimensions;
% directions are the columns and follows the DICOM convention;
%
% if dfile is empty, the k data will not be transformed.
% ds the difference in fov centers between T2 and Fatnav (T2-fatnav)
% k2 will be in the ro pe par coord system.

Nc=size(Data,4);
sz=size(Data);
nro=size(Data,1);

Data=reshape(Data,[nro,sz(2)*sz(3),Nc]);
iRep=iRep(:);

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


xform=afni_motionPar2Mat(dfile);
if max(iRep(:))~=size(xform,1)
    error('max(iRep(:))~=size(xform,1)');
end

xform=reshape(xform,[size(xform,1),4,3]);
xform=permute(xform,[3,2,1]);


k2=0*k;

for j=1:size(Data,2)
    
    ind2= (j-1)*nro+1:j*nro;
    ix=iRep(j);
    T= TransformMatrix_2CoordSystems(eye(3),ro_pe_par); % xform matrix in the ro-pe-par coordinator system
    xform_new=TransformMatrix_NewCoordSystem(xform(:,1:3,ix),T);
    xform_new(:,4)=T*xform(:,4,ix);
    
    [k2(ind2,:),Data(:,j,:)] = kspace_xform(k(ind2,:),Data(:,j,:),xform_new,voxsize,T*ds(:));
end

k2(:,2)=k2(:,2)*voxsize(2)/voxsize(1); %scale back to the original voxel size; 12/31/2019
k2(:,3)=k2(:,3)*voxsize(3)/voxsize(1); %scale back to the original voxel size

%%
Data = reshape(Data,[length(Data(:))/Nc,Nc]);


function rotmat=get_ro_pe_par(a)


norm=get_norm(a);

dphi=0; %don't know how to get the inplane rotation angle
rotmat=NormInplaneRot2Rotmat(norm,dphi);
rotmat(:,3)=-rotmat(:,3);
rotmat(:,2)=-rotmat(:,2);


function [kscaled,d3] = kspace_xform(k,data,xform,voxsize,ds)
% phase encoding is the 2nd dimension
% k and kscaled: (nro*npe)*3; -0.5 0.5; in the ro-pe-par system
% data: nro*npe*nch
% xform: 3*4; in the ro-pe-par system;
% ds: 1*3 or 3*1; in the ro-pe-par system
% voxel size along [ro, pe, and par]; 1*3

%k2=xform(:,1:3)*k';
%k2=k2';

k2=xform(:,1:3)*(k./voxsize*2*pi)';

k2=k2';

kscaled=k2.*voxsize/2/pi;

k3=reshape(k2',[3,size(data,1),size(data,2)]);
k3=repmat(k3,[1,1,1,size(data,3)]);

k3=permute(k3,[2,3,4,1]);
shift=xform(:,4)'; % phase shift
%shift=shift-((xform(:,1:3)-eye(3))*ds(:))';  % result is worse with - sign
shift=shift+((xform(:,1:3)-eye(3))*ds(:))';

shift=reshape(shift,[1,1,1,3]);
shift=repmat(shift,[size(data),1]);
d3 = data.*exp(-1i*sum(shift.*k3,4));

function norm=get_norm(a)
   dSag=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dSag');
   dCor=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dCor');
   dTra=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dTra');
 
    norm=[dSag,dCor,dTra];



function pos=get_pos(a)


 if isfield(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1},'sPosition')
     pos(1)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dSag');
     pos(2)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dCor');
     pos(3)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dTra');   
   else
     pos=[0,0,0];
   end
 
function res=get_field(s,fname)
    if isfield(s,fname)
      res=getfield(s,fname);
   else
    res=0;
   end



