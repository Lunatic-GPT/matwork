function Calc_MagMoment_Neurite(ph,mag,V,TE,B0,vox_size,prefix,vox_size_new,interp)

%B0: Telsa
%vox_size: 1*3 in mm
%vox_size_new: in mm;
% TE: ms
% V: 4 dim where the last dim has length 3 giving the direction of the DMV
% (in i,j,k coord system)
% at that voxel


if ~exist('vox_size_new','var')

vox_size_new = 0.2;
end

if ~exist('interp','var')
  %interp=ceil(vox_size(1:3)*2/vox_size_new);
  
  interp=2;
end
myAffine = true;

%ph=ri_mat_dcm(ph,[]);
%magdcm=ri_mat_dcm(mag,[]);
ph=double(ri(ph));
magdcm=double(ri(mag));

%%
B0_dir=[0,0,1];

phrange=max(ph(:))-min(ph(:));

data=double(magdcm).*exp(1i*double(ph)/double(phrange)*2*pi-1i*pi);

np=5;  %should be odd;


dchi =0.45; % ppm, SI unit

rad_mom=[0.4         1.2];
i_rad_mom=round(rad_mom/vox_size_new);

mask_calc=any(V>0,4);
%mask_calc=mask2>0&fr>=0.1*max(fr(mask2>0));

npix=sum(vec(mask_calc));

dim=round(25*0.4/vox_size_new);
mag_all=zeros(dim,dim,dim,npix,'single');
ph_all=zeros(dim,dim,dim,npix,'single');
deg=zeros(1,npix);

dchi_a2=zeros(np,np,np,npix,'single');
resid=zeros(np,np,np,npix,'single');

count=0;

for i=1:size(V,1)
    for j=1:size(V,2)
        for k=1:size(V,3)
            tic;
            if mask_calc(i,j,k)==0
                continue;
            end
   
            count=count+1;
            sub=[i,j,k];
            vd=squeeze(V(i,j,k,:));
         %   vd(3)=-vd(3);
            p1=sub+vd'*0.5;
            p2=sub-vd'*0.5;
            tic;
            [mag_tmp,ph_tmp,deg(count)]=rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,vox_size(1:3),data,32*ones(1,3),interp,vox_size_new,[],myAffine);
            
            disp(toc);
            mag_all(:,:,:,count)=crop(mag_tmp,[dim,dim,dim]);
            ph_all(:,:,:,count)=crop(ph_tmp,[dim,dim,dim]);
            
            for i2=1:np
                for j2=1:np
                    
                    for k2=1  %do not over count
                        center=ceil((size(ph_all(:,:,:,1))+1)/2)+[i2,j2,k2]-(floor(np/2)+1)*[1,1,1];
                        m_momtmp = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(1),center,1);
                        m_momtmp2 = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(2),center,1);
                        mask = m_momtmp==0&m_momtmp2>0;
                     
                        [dchi_a2(i2,j2,k2,count),resid(i2,j2,k2,count),cc(i2,j2,k2,count)]=magMoment_outSidePattern(ph_all(:,:,center(3),count)*180/pi,TE,B0,deg(count),mask,center(1:2),vox_size_new*[1,1]);  %180/pi added 11/3/2017
                        
                        dchi_a2(i2,j2,k2,count)=-dchi_a2(i2,j2,k2,count); % since the first dimension is along B0 instead of the second dim assumed in magMoment_outSidePa
%                         if debug
%                             disp(sqrt(dchi_a2(i2,j2,k2,count)/dchi)*1e3);
%                         end
                    end
                end
            end
            
            fprintf('Time remaining = %f s\n',toc*(npix-count));
        end
    end
end


dchi_a2(dchi_a2<0)=0;

dchi_a2=reshape(dchi_a2,[np*np*np,npix]);

y=zeros(1,npix);

ycc=zeros(1,npix);
imax_cc=zeros(1,npix);

for i=1:npix
   [ycc(i),indm]=max(vec(cc(:,:,:,i)));
    y(i) = dchi_a2(indm,i);
    imax_cc(i)=indm;
end

ccimg=single(0*mask_calc);
ccimg(mask_calc>0)=ycc;

rad = single(0*mask_calc);
rad(mask_calc>0)=sqrt(y/dchi)*1000;  % in units of um

%save(prefix,'dchi_a2','deg','dchi','rad','mask_calc','ccimg','cc','np','imax_cc','mag_all','ph_all');
 save(prefix,'dchi_a2','deg','dchi','rad','mask_calc','ccimg','cc','np','imax_cc');
   






