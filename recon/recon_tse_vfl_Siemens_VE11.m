function kdata=recon_tse_vfl_Siemens_VE11(fname,prefix_out)
% fname: *.mat file
% should contain Data, Line, Partition
% use the Siemens grappa setting, i.e., do IFT on both fully sampled
% dimensions before GRAPPA.
% 9/28/2020: add fpro parameter, was using prefix/prefix.pro
% suitable for data with transvere and inplane rotation = 0;

a=mapVBVD(fname);

lMatrix(1)=a{2}.hdr.Config.NImageCols;
lMatrix(2)=a{2}.hdr.Config.NImageLins;
lMatrix(3)=a{2}.hdr.Config.NoImagesPerSlab;
lMatrix(4)=a{2}.hdr.Config.NoOfFourierPartitions;
imSize=lMatrix([1,2,4]);

%imSize(1)=a{2}.hdr.Config.NImageCols;
%imSize(2)=a{2}.hdr.Config.NImageLins;
%imSize(3)=a{2}.hdr.Config.NImagePar;

d=readMeasDat_savemat(fname);

[Data3,Slines]=sort_kData(d.Data,d.Line,d.Partition,imSize);    

for i=1:size(Data3,4)
   Data3(:,:,:,i)=ifft1c(Data3(:,:,:,i),3);
end

[res,kdata]=recon_grappa2D_Siemens(Data3,Slines,imSize(2));

%coef_ch0=reshape(squeeze(coef(:,:,:,1,:)),[size(coef,1),Nc,size(coef,2)/Nc,size(coef,3),size(coef,5)]);
%save(fname, 'res','coef','-v7.3');% coef takes a lot of memory
%save([prefix_out,'.mat'], 'res','-v7.3');

%save tmp res kdata
o=get_rotmat(a,size(res));
mat2niigz(res,'',prefix_out,false,o);

if nargout==0
    kdata=[];
end

%
%    
% 
% for i=1:size(kdata,4)
%    fname=sprintf('%s_GRAPPA_Siem_Ch%d.mat',prefix,i); 
%    kdata_sc=kdata(:,:,:,i);
%    save(fname,'kdata_sc');  
% end
function o=get_rotmat(a,sz)
     fov(1)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
     fov(2)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
     fov(3)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
     
   if isfield(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1},'sPosition')
     pos(1)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dSag');
     pos(2)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dCor');
     pos(3)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dTra');   
   else
     pos=[0,0,0];
   end
     o.voxsize=fov./sz(:)';
     o.center=pos;
     
     o.rotmat=get_ro_pe_par(a);     
     o.orient='SPR';  
     o.sz=sz;
     o=center2pos_o(o,o.center);

function rotmat=get_ro_pe_par(a)

   dSag=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dSag');
   dCor=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dCor');
   dTra=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dTra');
 
    norm=[dSag,dCor,dTra];


dphi=0; %don't know how to get the inplane rotation angle
rotmat=NormInplaneRot2Rotmat(norm,dphi);
rotmat(:,3)=-rotmat(:,3);
rotmat(:,2)=-rotmat(:,2);


function [Data2,Slines]=sort_kData(Data,Line,Partition,imSize)


%p=unique(Line);

p=Line; 
s=Partition;

Nc=size(Data,2);

p_unq=unique(p);

Data2 = zeros([imSize(1),length(p_unq),imSize(3),Nc],'single');

for i=1:length(s)
    Data2(:,p(i)==p_unq,s(i),:)=Data(:,:,i);
end

Slines=p_unq;


function res=get_field(s,fname)
    if isfield(s,fname)
      res=getfield(s,fname);
   else
    res=0;
   end
