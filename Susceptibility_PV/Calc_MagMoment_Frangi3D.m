function Calc_MagMoment_Frangi3D(ph,mag,mask,TE,B0,vox_size)

%ph='TE2';
%mag='../MAG_IMAGES_0033/TE2';
%mask='mask_wm.nii';


%mask='mask_wm.nii';
b=load_nii(mask);
mask=b.img;
mask2=permute(mask,[2,1,3]);


%%

ph=ri(ph,[]);
magdcm=ri(mag,[]);
%%
options.FrangiScaleRange= [0.2 2];
options.FrangiScaleRatio = 0.2;
options.BlackWhite=false;
fname=['frangi_',num2str(options.FrangiScaleRange(1)),'_',num2str(options.FrangiScaleRange(2)),'_',num2str(options.FrangiScaleRatio),'.mat'];
if ~exist(fname,'file')
    [fr,whatScale,Voutx,Vouty,Voutz]=mFrangiFilter3D(ph,options,'',mask2);
    
    V=cat(4,Voutx,Vouty,Voutz);
    phi=atan2(V(:,:,:,2),V(:,:,:,1))*180/pi;
    th=acos(V(:,:,:,3)./sqrt(sos(V,4)))*180/pi;
    fr=single(fr);
    phi=single(phi);
    th=single(th);
    V = single(V);
    
    save(fname, 'fr','V','phi','th');
    
else
    %%
    fr=ri('frangi.mat');
    V=ri('frangi_vesselDir.mat');
end
%%

%clusterize2(fr>0.1,6);
B0_dir=[0,0,1];

vox_size_new=0.2;

phrange=max(ph(:))-min(ph(:));

interp=ceil(vox_size(1:3)*2/vox_size_new);
data=double(magdcm).*exp(1i*double(ph)/double(phrange)*2*pi-1i*pi);

np=1;  %should be odd;


dchi =0.45; % ppm, SI unit

rad_mom=[0.4         1.2];
i_rad_mom=round(rad_mom/vox_size_new);

mask_calc=mask2>0&fr>=0.1*max(fr(mask2>0));

npix=sum(vec(mask_calc));

mag_all=zeros(50,50,50,npix,'single');
ph_all=zeros(50,50,50,npix,'single');
deg=zeros(1,npix);

dchi_a2=zeros(np,np,np,npix,'single');
resid=zeros(np,np,np,npix,'single');

count=0;

for i=1:size(fr,1)
    for j=1:size(fr,2)
        for k=1:size(fr,3)
            tic;
            if mask2(i,j,k)==0 || fr(i,j,k)<0.1*max(fr(mask2>0))
                continue;
            end
            
            count=count+1;
            sub=[i,j,k];
            vd=squeeze(V(i,j,k,:));
            vd(3)=-vd(3);
            p1=sub+vd'*0.5;
            p2=sub-vd'*0.5;
            
            [mag_tmp,ph_tmp,deg(count)]=rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,vox_size(1:3),data,100*ones(1,3),interp,vox_size_new);
            mag_all(:,:,:,count)=crop(mag_tmp,[50,50,50]);
            ph_all(:,:,:,count)=crop(ph_tmp,[50,50,50]);
            
            for i2=1:np
                for j2=1:np
                    
                    for k2=1:np
                        center=ceil((size(ph_all(:,:,:,1))+1)/2)+[i2,j2,k2]-(floor(np/2)+1)*[1,1,1];
                        m_momtmp = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(1),center,1);
                        m_momtmp2 = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(2),center,1);
                        mask = m_momtmp==0&m_momtmp2>0;
                        [dchi_a2(i2,j2,k2,count),resid(i2,j2,k2,count)]=magMoment_outSidePattern(ph_all(:,:,center(3),count),TE,B0,deg(count),mask,center(1:2),vox_size_new*[1,1]);
                        dchi_a2(i2,j2,k2,count)=-dchi_a2(i2,j2,k2,count); % since the first dimension is along B0 instead of the second dim assumed in magMoment_outSidePa
                        if debug
                            disp(sqrt(dchi_a2(i,j,k,ipix)/dchi)*1e3);
                        end
                    end
                end
            end
            
        end
    end
end

%%
dchi_a2(dchi_a2<0)=0;

y=dchi_a2(1,1,1,:);
yres=resid(1,1,1,:);

rad = 0*mask_calc;
rad(mask_calc>0)=sqrt(y/dchi)*1000;


img_resid=single(0*mask_calc);
img_resid(mask_calc>0)=yres;


mask_calc=uint8(mask_calc);

save res_all_center mag_all ph_all deg  dchi rad img_resid mask_calc

% %% just for debug
% load res_all.mat
% 
% npix=sum(mask_calc(:)>0);
% %mag_all=mag_all2(:,:,:,2:end);
% %ph_all=ph_all2(:,:,:,2:end);
% vox_size_new = 0.2;
% 
% np=1;  %should be odd;
% 
% B0 = 7;
% TE = 15;
% dchi =0.45; % ppm, SI unit
% 
% dchi_a2=zeros(np,np,np,npix,'single');
% 
% resid=zeros(np,np,np,npix,'single');
% rad_mom=[0.4         1.2];
% i_rad_mom=round(rad_mom/vox_size_new);
% 
% debug=true;
% if debug
%     ind=sub2ind(size(mask_calc),216,236,70);
%     
%     ipix=find(find(mask_calc>0)==ind);
% end
% 
% 
% for ipix=1:npix
%     tic;
%     for k=1:np
%         for i=1:np
%             for j=1:np
%                 
%                 center=ceil((size(ph_all(:,:,:,1))+1)/2)+[i,j,k]-(floor(np/2)+1)*[1,1,1];
%                 m_momtmp = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(1),center,1);
%                 m_momtmp2 = mask_circle(size(ph_all(:,:,1,1)),i_rad_mom(2),center,1);
%                 mask = m_momtmp==0&m_momtmp2>0;
%                 [dchi_a2(i,j,k,ipix),resid(i,j,k,ipix)]=magMoment_outSidePattern(ph_all(:,:,center(3),ipix),TE,B0,deg(ipix),mask,center(1:2),vox_size_new*[1,1]);
%                 dchi_a2(i,j,k,ipix)=-dchi_a2(i,j,k,ipix); % since the first dimension is along B0 instead of the second dim assumed in magMoment_outSidePa
%                 if debug
%                     disp(sqrt(dchi_a2(i,j,k,ipix)/dchi)*1e3);
%                 end
%                 
%             end
%         end
%     end
%     toc;
% end
% 
% % %%
% % dchi_a2=reshape(dchi_a2,[np^3,575]);
% %
% % y=dchi_a2(14,:);
% % yres=resid(2,2,2,:);
% 
% dchi_a2(dchi_a2<0)=0;
% 
% y=dchi_a2(1,1,1,:);
% yres=resid(1,1,1,:);
% 
% rad = 0*mask_calc;
% rad(mask_calc>0)=sqrt(y/dchi)*1000;
% 
% 
% img_resid=single(0*mask_calc);
% img_resid(mask_calc>0)=yres;
% 
% 
% mask_calc=uint8(mask_calc);



%%













