%% to do: with 1TE, the vein intensity is negative; non-physical, needs to add constraint to the problem.
%% no longer used, use gui_spv instead.  When nullVein, results not stable.
swiname='swi_TE1';
dir_trace='NeuriteTracer';
prefix_trace = 'SWI_Images_Inv-exported';  %or the prefix of the vein mask.
maskfile='SWI_veinmask';
phdir='Pha_Images';
mgdir='Mag_Images';
vox_size_new=0.2; %unit mm.

B0_dir=[0,0,1];
B0=7; % field strength
crop_size=[64,64,32];  % crop size before interpolation.
lout = 32;   % output matrix size after rotate data.
lseg=3; % units mm.  path segment length.

step =1;
% step 1:  generate inverted nii file from nii saved in imageJ
% step 2: generate mask based on the neurite traces. check mask and save it as
% maskfile.mat; especially check for neighboring paths.
% step 3: determined paths from the mask; 
% step 4: rotate data
% step 5: calculate susceptibility
% step 6: disp result
  method = 'own';
%% step 1: load SWI dicom files into imageJ.   Save as analyze file

if any(step==1)
    d=load_nii(swiname);
    d.img=max(d.img(:))-d.img;
    save_nii(d,[swiname,'_Inv']);
    
    % Step 2: load _Inv.nii; draw traces in imageJ
end
    %% step 2
if any(step ==2)
    dcm2mat_siemens(swiname);
    add_neurite_traces([swiname,'.mat'],dir_trace,prefix_trace,{'xy'},-2);
    
    %check results in img
end    
    %% step 3
if any(step==3)
    
        ThinningPathFind([maskfile,'.mat']);
    
    %% step 3b
    load(['ThPth_',maskfile]);  % need to fix connected paths afterwards
    %% check path length
    for i=1:length(i_ind_path)
        len=pathLength(i_ind_path{i},ind{i},voxsize,size(c));
        fprintf('Path Length %d: %f\n',i,len);
    end
    
    
     pathCurved_nofix(['ThPth_',maskfile]);  %fixedpath

 
     fixGap(['ThPth_',maskfile,'_Curv']);   %probv

end 
    %% step 4: rotate data
if any(step==4)
    
    [p1,p2,an,resid_pathDir]=Path_Orientation(['ThPth_',maskfile,'_Curv_Dv.mat'],lseg);
    phdcm=ri(phdir,[]);
    magdcm=ri(mgdir,[]);
    phrange=max(phdcm(:))-min(phdcm(:));
    vox_size=ri(['ThPth_',maskfile,'_Curv_Dv.mat'],[],[],'voxsize');
    
    
    interp=ceil(vox_size(1:3)*2/vox_size_new);
    data=double(magdcm).*exp(1i*double(phdcm)/double(phrange)*2*pi-1i*pi);
    
    lmax=maxlen_cell(an);
    
    mag=zeros(lout,lout,lmax,length(an));
    
    ph=zeros(lout,lout,lmax,length(an));
    deg=cell(1,length(an));
    for i=3%1:length(an)
        for j=3%1:length(an{i})
            %  prefix=sprintf('rotateDataOutput/path_%d_seg_%d',i,j);
            [mag_tmp,ph_tmp,deg{i}(j)]=rotate_data_Line2dim3_B2dim1(p1{i}(j,:),p2{i}(j,:),B0_dir,vox_size(1:3),data,crop_size,interp,vox_size_new);
            c0=ceil(size(mag_tmp)/2+0.5);
            
            %c=centerFinder(mag_tmp(:,:,c0(3)).*exp(1i*ph_tmp(:,:,c0(3))),2,4);
            c=centerFinder(mag_tmp(:,:,c0(3)),1,4,1);
            
            
            mag(:,:,j,i) = mag_tmp(c(1)-lout/2:c(1)+lout/2-1,c(2)-lout/2:c(2)+lout/2-1,c0(3));
            ph(:,:,j,i) = ph_tmp(c(1)-lout/2:c(1)+lout/2-1,c(2)-lout/2:c(2)+lout/2-1,c0(3)); 
            
          %  mag(:,:,j,i) = mag_tmp(c0(1)-lout/2:c0(1)+lout/2-1,c0(2)-lout/2:c0(2)+lout/2-1,c0(3));
          %  ph(:,:,j,i) = ph_tmp(c0(1)-lout/2:c0(1)+lout/2-1,c0(2)-lout/2:c0(2)+lout/2-1,c0(3)); 
            
        end
    end
    
    ph_prefix=strtok(phdir);
    save(['Rot_Ph_Mag_',maskfile,'.mat'],'mag','ph','deg','resid_pathDir');
    
end
if any(step==5)
    %% step 5: calc susceptibility
    ph_prefix=strtok(phdir);
    
    TE=readdPar(phdir,'EchoTime');    
    
    
    mag=ri(['Rot_Ph_Mag_',maskfile,'.mat'],'','','mag');
    ph=ri(['Rot_Ph_Mag_',maskfile,'.mat'],'','','ph');
    cd = mag.*exp(1i*ph);
    
    deg=ri(['Rot_Ph_Mag_',maskfile,'.mat'],'','','deg');
    resid_pathDir=ri(['Rot_Ph_Mag_',maskfile,'.mat'],'','','resid_pathDir');
    
    res=cell(1,length(deg));
    flag=cell(1,length(deg));
    resid_magMoment=cell(1,length(deg));
    m_large=mag*0;
    m_vessel=mag*0;
    
    npath=length(deg);
  
    
    for i=1:npath
        for j=1:length(deg{i})
            disp([i,j]);
            c=ceil((size(mag)+1)/2);
            m_vessel(:,:,j,i)=mask_circle(size(mag),2,c(1:2),1);
            if strcmp(method,'mri')
                tm1=mask_circle(size(mag),4,c(1:2),1);
                tm2=mask_circle(size(mag),8,c(1:2),1);
                m_large(:,:,j,i) = tm1;% tm2>0&tm1==0;
                [res{i}(j,:),flag{i}(j),resid_magMoment{i}(j)]= calc_vessel_susceptibility(cd(:,:,j,i),m_vessel(:,:,j,i),m_large(:,:,j,i),deg{i}(j),TE,B0,vox_size_new);
            else
                [res{i}(j,:),flag{i}(j),resid_magMoment{i}(j)]= calc_vessel_susceptibility_nullVein(cd(:,:,j,i),m_vessel(:,:,j,i),vox_size_new*(5:15),deg{i}(j),TE,B0,vox_size_new);
                
                
            end
        end
    end
    if strcmp(method,'mri')
      mname=sprintf('masks_vessels_%s_mri.mat',maskfile);
      save(mname,'m_vessel','m_large');
      resname=sprintf('S_PV_%s_mri.mat',phdir);
     
    else
      mname=sprintf('masks_vessels_%s_own.mat',maskfile);
      save(mname,'m_vessel');
      resname=sprintf('S_PV_%s_own.mat',phdir);
     
    end
     save(resname,'res','vox_size_new','TE','deg','B0','flag','resid_magMoment','resid_pathDir');
end
if any(step==6)
    if strcmp(method,'mri')
        resname=sprintf('S_PV_%s_mri.mat',phdir);
    else
        resname=sprintf('S_PV_%s_own.mat',phdir);
    end
    load(resname);
    
    tmp=cell2array(res,1);
    chi=tmp(:,1);
    r=tmp(:,2);
    mom=tmp(:,3);
    
    flga=cell2array(flag);
    rsdp=cell2array(resid_pathDir);
    rsdm=cell2array(resid_magMoment);
    
    figure;plot(chi(flga==1),r(flga==1),'ro');
    hold on;
    plot(chi(flga~=1),r(flga~=1),'k>');
    
    figure;plot(chi(flga==1),r(flga==1),'ro');
    hold on;
    plot(chi(flga~=1),r(flga~=1),'k>');
    
    thr=0.06;
    figure;plot(chi(rsdp>thr),r(rsdp>thr),'ro');
    hold on;
    plot(chi(rsdp<=thr),r(rsdp<=thr),'k>');
    
    thr=0.1;
    figure;plot(chi(rsdm>thr),r(rsdm>thr),'ro');
    hold on;
    plot(chi(rsdm<=thr),r(rsdm<=thr),'k>');
    
%%
figure; 
subplot(1,2,1);

xc=0.3:0.02:0.6;
%hist(chi(tmp(:,5)>0),xc);
hist(chi(chi<1),xc);
xlim([min(xc),max(xc)]);
set(gca,'FontSize',14);
xlabel('Susceptibility (ppm)');
ylabel('Number of Segments');

subplot(1,2,2);
xc=0:0.01:0.4;
hist(real(r),xc);
xlim([min(xc),max(xc)]);
set(gca,'FontSize',14);
xlabel('Radius (mm)');
ylabel('Number of Segments');

    
end   



%% note: path 3; seg 3, negative chi when m_large radius = 6, because phase and mag centers do not match?





