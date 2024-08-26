function recon_ESPIRiT3D_FatNav(fname,dfile,njobs,nmapjobs,file_res0,bigmem)

% fname='meas_MID85_tse_vfl_pss_FatNav_FID23810.mat';
% dfile = 'Motion_MID85.1D';
% fov = [210.0,165.7,97.6];  % FOV before accounting for the oversampling
% factor
% lMatrix = [512,405,244,256];
% 8/25/2019: remove fov and lMatrix parameters; 

while ~exist(dfile,'file') || ~exist(fname,'file')
    pause(600);
end


load(fname);


prefix=strtok2(fname,'.');
    fov(1)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dReadoutFOV');
    fov(3)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dThickness');
    fov(2)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dPhaseFOV');
    
    
    lBaseResolution = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lBaseResolution');
    lPhaseEncodingLines =readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lPhaseEncodingLines');
    lPartitions = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lPartitions');
    lImagesPerSlab = readsPar(fullfile(newdname,prefix,[prefix,'.pro']),'lImagesPerSlab');
    
lMatrix=[lBaseResolution,lPhaseEncodingLines,lImagesPerSlab,lPartitions];

        pnav_sag = readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dSag');
        pnav_cor= readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dCor');
        pnav_tra=readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dTra');      
        if isempty(pnav_sag)
            pnav_sag=0;
        end
        if isempty(pnav_cor)
            pnav_cor=0;
        end
        if isempty(pnav_tra)
            pnav_tra=0;
        end
           
        p_sag = readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dSag');
        p_cor= readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dCor');
        p_tra=readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dTra');
        
        
         if isempty(p_sag)
            p_sag=0;
        end
        if isempty(p_cor)
            p_cor=0;
        end
        if isempty(p_tra)
            p_tra=0;
        end
        
Data = fft1c(Data,1);   % to k space;
nro=size(Data,1);
Nc=32;
% 
disp('calculate maps');
fov_os=[fov(1:2),fov(3)*lMatrix(4)/lMatrix(3)];  % fov_os after oversampling along Partition direction
ro_pe_par=[1,-2,-3];
ds = [p_sag,p_cor,p_tra]-[pnav_sag,pnav_cor,pnav_tra];
[k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov_os,lMatrix([2,4]),ro_pe_par,ds);


%prep_bsub_ESPIRiT3D_maps(fname,dfile,fov_os,lMatrix([2,4]),nmapjobs,ncomp);
%% old prep_bsub_ESPIRiT3D_maps
    nmaps = 2;
    lRefLinesPE=24; % fix for now
    kthr=54/lMatrix(2)/2;
   
    prefix2 = [prefix,'_kCalib_maps'];
    ind=linspace(0,nro,nmapjobs+1);
    ind=round(ind);
    ind1=ind(1:end-1)+1;
    ind2=ind(2:end);
    
    if wait_for_files(prefix2,ind1,ind2,0,0)  %if not found
        prep_bsub_ESPIRiT3D_maps(k2,Data,Nc,lMatrix,kthr,lRefLinesPE,prefix,nmapjobs,nmaps);
        if ispc
            wait_for_files(prefix2,ind1,ind2,0);
        else
            wait_for_files(prefix2,ind1,ind2);
        end
    end
%%

if ~bigmem &&~ispc
    mid=strtok_no(prefix,'_',2);
    script_name=sprintf('recon_script_%s_2.sh',mid);
    fid=fopen(script_name,'w');
     fprintf(fid,'#!/bin/bash\n');
        
        fprintf(fid,'#SBATCH --job-name=recon\n');
        fprintf(fid,'#SBATCH --ntasks=1\n');
        fprintf(fid,'#SBATCH --time=24:00:00\n');
        fprintf(fid,'#SBATCH --mem=400000\n');
        fprintf(fid,'#SBATCH --partition=bigmem\n');
        fprintf(fid,'#SBATCH --qos bigmem_access\n');
               
        fprintf(fid,'matbgk "recon_ESPIRiT3D_FatNav(''%s'',''%s'',[%4.1f,%4.1f,%4.1f],[%d,%d,%d,%d],%d,%d,''%s'',true)" logrecon_%s\n',fname,dfile,fov,lMatrix,njobs,nmapjobs,file_res0,mid);
        fclose(fid);
        cmd=sprintf('sbatch %s',script_name);
        unix(cmd);
        
    return;
end
tmp=load([prefix,'_kCalib.mat'],'imSize');
imSize=tmp.imSize;



    fmap_ind = [prefix2,'_Ind.mat'];

    save(fmap_ind,'ind1','ind2','nro','Nc','nmaps');
    
    
 distMem = false;%njobs is not used;
 nufft=NUFFT3DDistMem(k2,1,[0,0,0],imSize,njobs,prefix2,fmap_ind,distMem);


prefix_save=unique_name([prefix,'_recon_MotionCorr']);

if ~exist('file_res0','var') || isempty(file_res0)
    
    disp('Calculating initial image');
    res=nufft'*Data;
    disp('Initial image calculation done');
    
    save([prefix_save,'_iter0.mat'],'res', '-v7.3');
    res0=res;
    
else
    res0=load(file_res0);
    res0=res0.res;
end



XOP = Wavelet('Daubechies',4,6);
disp('Start CGL1ESPIRiT');
  tcg=tic;
nIterCG = 5;

% 
% if combine_esp_nufft
  res = cgL1ESPIRiT_nomap(Data, res0, nufft, nIterCG,XOP,0,0,1);
% else
%   res = cgL1ESPIRiT(double(DATA), zeros(sx,sy,2), FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);
% end

disp('CGL1ESPIRiT finished');
toc(tcg);

if exist('file_res0','var') && ~isempty(file_res0)
    prefix_save=file_res0(1:end-4);
end
    save([prefix_save,'_iter',num2str(nIterCG),'.mat'], 'res','-v7.3');

   function prep_bsub_ESPIRiT3D_maps(k2,Data,Nc,lMatrix,kthr,lRefLinesPE,prefix,njobs,nmaps)
        
mid=strtok_no(prefix,'_',2);
nro=lMatrix(1);
        %%
        sel=abs(k2(:,2))<=kthr & abs(k2(:,3))<=kthr;
  k2_center=k2(sel,:);
  
imSize=lMatrix([1,2,4]);
imSize_lowres = round(imSize*kthr/0.5);%[length(pe_center),length(ro_center),length(par_center)];
imSize_lowres(1)=imSize(1);  % no change in the readout direction.
  
        Data = reshape(Data,[length(Data(:))/Nc,Nc]);
       
        Data=Data(sel,:);
        nufft=NUFFT3D(k2_center.*repmat(imSize./imSize_lowres,[size(k2_center,1),1]),1,[0,0,0],imSize_lowres,1,1);
        %   tic;
        im = nufft'*Data;
        
        
        kData=fft1c(fft1c(im,2),3);   % this matches with d2
        
        acs_par=size(kData,3)/2-11:size(kData,3)/2+12;
        
        acs_line=size(kData,2)/2-lRefLinesPE/2+1:size(kData,2)/2+lRefLinesPE/2;
        
        kCalib=squeeze(kData(:,acs_line,acs_par,:));
        
        
        fname2 = [prefix,'_kCalib.mat'];
        save(fname2, 'kCalib','imSize');
        
       
        
        ind=linspace(0,nro,njobs+1);
        
        ind=round(ind);
        
        for i=1:length(ind)-1
            
            fmapname=sprintf('%s_maps_%d_%d.mat',strtok(fname2,'.'),ind(i)+1,ind(i+1));  %do not regenerate
            if exist(fmapname,'file')
                continue;
            end
            
            if ispc
                calc_ESPIRiT3D_maps(fname2,ind(i)+1:ind(i+1),nmaps);%d:%d,%d) logmaps_%d\n',fname2,ind(i)+1,ind(i+1),nmaps,ind(i)+1)
            else
                
                
                fname3=unique_name(sprintf('sbatch_ESPIRiT3D_maps_%s_%d.sh',mid,ind(i)+1));
                
                fid = fopen(fname3,'w');
                
                fprintf(fid,'#!/bin/bash\n');
                
                fprintf(fid,'#SBATCH --job-name=map%d\n',i);
                fprintf(fid,'#SBATCH --ntasks=1\n');
                fprintf(fid,'#SBATCH --time=24:00:00\n');
                fprintf(fid,'#SBATCH --mem=32000\n');
                % fprintf(fid,'#SBATCH --partition=bigmem\n');
                %  fprintf(fid,'#SBATCH --qos bigmem_access\n');
                
                %fprintf(fid,'bsub -M 20 -o logmaps_%d.%%J matbgk "calc_ESPIRiT3D_maps(''%s'',%d:%d)" logmaps_%d\n',ind(i)+1,fname2,ind(i)+1,ind(i+1),ind(i)+1);
           
                fprintf(fid,' matbgk "calc_ESPIRiT3D_maps(''%s'',%d:%d,%d)" logmaps_%s_%d\n',fname2,ind(i)+1,ind(i+1),nmaps,mid,ind(i)+1);
                
                fclose(fid);
                cmd=sprintf('sbatch %s',fname3);
                disp(['submit ',fname3]);
                pause(60);  % pause so that all scripts will be submitted properly, hopefully.
                unix(cmd);
            end
        end
        


function [k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov,lMatrix,ro_pe_par,ds)
% fov [3]: field of v for [ro, pe,par]
% lMatrix: matrix size for [pe,par];
%
% ro_pe_par: [3]; first,second, third elements for ro, pe, par, respectively,
                 % 1 - x (LR); 2 - y (AP); 3 - z (IS);
if ~exist('ro_pe_par','var')
    ro_pe_par=[2,1,3];
end


xform=afni_motionPar2Mat(dfile);  % should get the same value with the following command  
  
xform=reshape(xform,[size(xform,1),4,3]);
xform=permute(xform,[3,2,1]);

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

if size(xform,3)~=length(unique(s))
    error('number of transformed wrong');
end

ns=length(unique(s2));
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
shift=shift+((xform(:,1:3)-eye(3))*ds(:))';

shift=shift./voxsize*2*pi;

shift=reshape(shift,[1,1,1,3]);


shift=repmat(shift,[size(data),1]);


d3 = data.*exp(-1i*sum(shift.*k3,4));  % 


