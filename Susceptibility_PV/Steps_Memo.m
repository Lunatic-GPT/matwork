% steps:
% 1. load dicom files into ITK-Snap;   
% the following true for FRH01\Mag_Images
% save main image as .nii
% x,y,z in ITK-snap corresponds to 1st,2nd, and 3rd dimensions in .nii
% data.  While in dicom file, [x,y,z] <-> [2nd, 1st, 3rd dimensions].
ph=load_nii('Pha_Images.nii');
mag=load_nii('Mag_Images.nii');

data=double(mag.img).*exp(1i*double(ph.img)/8192*2*pi);

B0_dir=[0,0,1];

vox_size=ph.hdr.dime.pixdim(2:4);

p1=[135,302,143];
p2=[119,302,145];

crop_size=[64,64,32];
interp=4;

vox_size_new=0.25;
str=sprintf('_%d',p1);
prefix=['TE20ms',str];

%%
rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,vox_size,data,crop_size,interp,vox_size_new*ones(1,3));

%%
pos_str1=sprintf('%d_',p1);
pos_str2=sprintf('%d_',p2);
pos_str1=pos_str1(1:end-1);
pos_str2=pos_str2(1:end-1);


d=sprintf('affine_%s-%s.mat',pos_str1,pos_str2);
m_large=sprintf('affine_%s-%s_Largeroi.mat',pos_str1,pos_str2);
m_vessel=sprintf('affine_%s-%s_Vesselroi.mat',pos_str1,pos_str2);

TE=readdPar('Pha_Images', 'EchoTime');
theta_B0=ri(d,'','','deg');
[res,flag]= calc_vessel_susceptibility(d,m_vessel,m_large,theta_B0,TE,7,vox_size_new);





