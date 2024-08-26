function recon_fl_fq_inplanemoco(fsb,interp_factor,prefix,Tfile,ind_inc)

% No retro-gating.
% Can work for both short and long term averaged data and data with multiple measurements. 
% Tfile: mat file containing the transformation matrices T; if empty then,
%        no moco will be applied; T should be of length num of average.
%        The will be applied to the images after applying the
%        interp_factor.
% ind_inc: mat file containing indices (variable ind_inc) to include in the average, 1 based, Should be of the
%          same length as T.

if ~exist('interp_factor','var')|| isempty(interp_factor)
    interp_factor=[2,3.333];
end

prefix_fsb=strtok2(filename(fsb),'.');
pro=fullfile(strtok2(fsb,'.'),[prefix_fsb,'.pro']);

%[dsb,lin,par,sl,ushSet,timeStamp,freePara,dfn,lin_dfn,par_dfn,set_dfn,dref_dfn,d_ps] = readMeasDat(fsb,inf,0,true);
dsbStruct = readMeasDat_savemat(fsb,Inf);
ushSet=dsbStruct.Set;
lin=dsbStruct.Line;
dsb=dsbStruct.Data;


pat=readsPar(pro,'lAccelFactPE');
npe=readsPar(pro,'lPhaseEncodingLines');
%%

nvenc=max(ushSet)+1;


nphase=Inf;
physio=struct;
freePara=zeros(length(lin),1);

dsb2=reorder_fl_fq_data(dsb,lin,pro,nvenc,freePara,nphase,physio,false,90);

if ~exist('ind_inc','var')
    ind_inc=1:size(dsb2,6);
else
    tmp=load(ind_inc);
    ind_inc=tmp.ind_inc;
end

dsb2=dsb2(:,:,:,:,:,ind_inc);  %[nro,lPhaseEncodingLines,nsl,32,nvenc,nphase];

if pat>1
    dsb2_mean=mean(dsb2(:,:,:,:,1,:),6);
    lin2_unq=sort(unique(lin(:)));
    for i=1:size(dsb2,6)
        [~,dsb2(:,:,:,:,:,i)]=recon_grappa2D(fft1c(dsb2(:,lin2_unq+1,:,:,:,i),1),lin2_unq+1,npe,fft1c(dsb2_mean(:,lin2_unq+1,:,:),1),lin2_unq+1);
    end
else
    dsb2= fft1c(dsb2,1);
end




im_sb = do_recon(dsb2);

orient=proDimCenter(pro,round(size(im_sb(:,:,1,1)).*interp_factor));


im_sb2=do_interp_dim12_local(im_sb,interp_factor);


if exist('Tfile','var') && ~isempty(Tfile)
    tmp=load(Tfile);
    im_sb3=apply_moco(im_sb2,tmp.T(ind_inc));
else
    im_sb3=im_sb2;
end
save_nii(mean(im_sb3,6),prefix,dsb(:,1:32),orient);


function im=apply_moco(im0,T)

im=0*im0;
tic;
for i=1:length(T)
    for i4=1:size(im0,4)
        for i3=1:size(im0,3)
            for i5=1:size(im0,5)
                for i6=1:size(im0,6)
                    im(:,:,i3,i4,i5,i6)=imwarp(im0(:,:,i3,i4,i5,i6),T(i),'OutputView',imref2d(size(im0(:,:,1,1,1,1))));
                end
            end
        end
        toc;
        time_left(i4,size(im0,4),toc,i,length(T));
    end
end


function d2=do_interp_dim12_local(d,interp)

sz=size(d);

sz(1:2)=round(sz(1:2).*interp);


nt=sz(6);
d2=zeros(sz,'single');
if any(interp>1)
    for l=1:nt
        for j=1:sz(4)
            for k=1:sz(5)
                fd=single(fft2c(d(:,:,:,j,k,l)));
                tmp=zpad(fd,sz(1),sz(2),sz(3));
                d2(:,:,:,j,k,l)=ifft2c(tmp);
            end
        end
    end
else
    d2=d;
end

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


function save_nii(d,prefix,nois,orient)
%do_interp_dim12(d,interp,mid,nois,suffix,do_scale,roi)
% im should be nro*npe*nsl*nch*nvenc*nt

% nois: nro*nch; can be empty
% phase_Only: only calculate the phase images, no mag or PC.
% 9/27/2020:
% now save in nii format
% deleted phase_only parameter;
% delete mid; now use prefix;
% add orient parameter;


nt=size(d,6);
max_mag=zeros(1,nt);

sz=size(d);

mag=zeros(sz(1),sz(2),sz(3),sz(5),nt,'uint16');
pc=zeros(sz(1),sz(2),sz(3),sz(5)-1,nt,'int16');
do_scale=true;

for l=1:nt
    [mag(:,:,:,:,l),pc(:,:,:,:,l),max_mag(l)]=coil_combine_PC(d(:,:,:,:,:,l),nois,false,do_scale);
end

if do_scale
    for l=1:nt
        tmp=mag(:,:,:,:,l);
        mag(:,:,:,:,l)=uint16(single(tmp)*max_mag(l)/max(max_mag));
    end
end

nro=size(mag,1);
nsl=size(mag,3);
nvenc=size(mag,4);

mag=reshape(mag,[nro,size(mag,2),nsl,nvenc*nt]);
pc=reshape(pc,[nro,size(pc,2),nsl,nt]);
save_nii_local(mag,orient,['Mag_',prefix]);
save_nii_local(pc,orient,['Phase_',prefix]);


function save_nii_local(d,orient,prefix)

nii=make_nii(d);
nii=nii_from_orient(nii,orient);
save_untouch_niigz(nii,prefix);

