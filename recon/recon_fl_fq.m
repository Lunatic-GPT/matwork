function recon_fl_fq(fsb,interp_factor,nphase,prefix,physio,DownsampleFactor,do_scale)

% if nphase is Inf; then set it equal to number of averages and do not use
% the physio data for retro-gating.
%physio: can be a file name; 
% works for both short-term and long term averaged data, or no average but
% multiple measurements.

if ~exist('DownsampleFactor','var') || isempty(DownsampleFactor) %downsample by crop the k-space.
    DownsampleFactor=[1,1];
end

if ~exist('interp_factor','var')|| isempty(interp_factor)
    interp_factor=[2,3.333];
end

if ~exist('nphase','var')
    nphase=12;
end

if ~exist('physio','var')  || isempty(physio)
    physio=struct;
end

if isa(physio,'char')
    load(physio);
end

if isfield(physio,'maxRate')
    maxRate=physio.maxRate;
else
    maxRate=80;
end

if ~exist('do_scale','var')
    do_scale=true;
end

prefix_fsb=strtok(filename(fsb),'.');
prefix_dir=fullfile(prefix_fsb,prefix_fsb);
if ~exist([prefix_fsb,'.mat'],'file')
    %[dsb,lin,par,sl,ushSet,timeStamp,freePara,dfn,lin_dfn,par_dfn,set_dfn,dref_dfn,d_ps] = readMeasDat(fsb,inf,0,true); 
    [dsbStruct,~,freePara] = readMeasDat(fsb,Inf,0,true);
     ushSet=dsbStruct.Set;
        lin=dsbStruct.Line;
        dsb=dsbStruct.Data;
   % save([prefix,'.mat'],'dsb','lin','par','sl','ushSet','freePara','dfn','lin_dfn','par_dfn','set_dfn','dref_dfn','d_ps');
else   
    load([prefix_fsb,'.mat']);
    if exist('Set','var')  % new format
        ushSet=Set;
        lin=Line;
        dsb=Data;
        clear Data;
    end
end

if exist('no_recon','var') && no_recon
    return;
end

nave=readsPar([prefix_dir,'.pro'],'lAverages');
seg=readsPar([prefix_dir,'.pro'],'lSegments');

if isempty(nave)
    nave=1;
end

pat=readsPar([prefix_dir,'.pro'],'lAccelFactPE');
npe=readsPar([prefix_dir,'.pro'],'lPhaseEncodingLines');
%%

nch=nChan_SiemensProt([prefix_dir,'.pro']);
nvenc=max(ushSet)+1;

%%

try
    phaseStab=readsPar([prefix_dir,'.pro'],'ucPhaseStabilize');   
    if strcmp(phaseStab{1},'0x1')
        phaseStab=true;
    else
        phaseStab=false;
    end
catch
    phaseStab=false;
end

old_way=1; % old way has no interpolation
if old_way
  dsb2=reorder_fl_fq_data(dsb,lin,prefix_dir,nvenc,freePara(:,4),nphase,physio,false,maxRate);
else
  ns=1;
  istart=find(freePara(1:nch:end,2)>0);
  dsb=reshape(dsb,[size(dsb,1),nch,ns,size(dsb,2)/nch/ns]);
  dsb=permute(dsb,[1,4,3,2]);
  dsb2=reorder_fl_fq_data_interp(dsb,istart,Set(1:ns*nch:end),nphase);
end
  
 
%%

if pat>1
    dsb2_mean=mean(dsb2(:,:,:,:,1,:),6);
    lin2_unq=sort(unique(lin(:)));
    for i=1:size(dsb2,6)
        %dsb2(:,:,:,:,:,i)
        [im,dsb2(:,:,:,:,:,i)]=recon_grappa2D(fft1c(dsb2(:,lin2_unq+1,:,:,:,i),1),lin2_unq+1,npe,fft1c(dsb2_mean(:,lin2_unq+1,:,:),1),lin2_unq+1); 
    end
else  
    dsb2= fft1c(dsb2,1); 
end


sz=size(dsb2);
sz(1:2)=round(sz(1:2).*DownsampleFactor);
dsb2=crop(dsb2,sz);

im_sb = do_recon(dsb2);
%[tmp,mid]=strtok(prefix,'_');
%mid=strtok(mid(2:end),'_');
pro=[prefix_dir,'.pro'];
orient=proDimCenter(pro,round(size(im_sb(:,:,1,1)).*interp_factor));

do_interp_dim12(im_sb,interp_factor,prefix,dsb(:,1:32),orient,do_scale);

%% add
% try
%   
%     for i=1:length(outname)
%         mat_addv(fullfile(outname(i).folder,outname(i).name),'voxsize',coord.voxsize);
%         mat_addv(fullfile(outname(i).folder,outname(i).name),'center',coord.center);
%         mat_addv(fullfile(outname(i).folder,outname(i).name),'rotmat',coord.rotmat);
%         mat_addv(fullfile(outname(i).folder,outname(i).name),'pos',coord.pos);
%         mat_addv(fullfile(outname(i).folder,outname(i).name),'orient',coord.orient);
%     end
%     
% catch
%     return;
% end


function im_sb = do_recon(dsb2)

tmp=sos(dsb2(:,:,:,:,1,1),4);
[tmp2,ind_max]=max(tmp(:));
ind_max=ind2subb(size(tmp),ind_max);
negphase=ind_max(1)-1;
%%
sz=size(dsb2);
im_sb=single(zeros(sz));

for i=1:size(dsb2,4)
    for j=1:size(dsb2,5)
        for k=1:size(dsb2,3)
            for iphase=1:size(dsb2,6)
                tmp=dsb2(:,:,k,i,j,iphase);
               
               if abs((sz(1)/2-negphase)/sz(1))>0.1 % no partial Fourier recon if less than 10%
                    tmp=cat(1,zeros(sz(1)/2-negphase,size(tmp,2)),tmp(1:negphase+sz(1)/2,:));
                    im_sb(:,:,k,i,j,iphase)=partialFT_pocs(tmp,40,true);     
                else 
                    im_sb(:,:,k,i,j,iphase)=ifft1c(ifft1c(tmp,1),2);
                end 
            end 
        end        
    end
end







  
  
