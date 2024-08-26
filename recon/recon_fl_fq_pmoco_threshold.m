function recon_fl_fq_pmoco_threshold(fname,interp_factor,prefix,irep,tran_threshold,rot_threshold)
%irep: the number of repetitions used for average calculation.
% tran_threshold in mm
% rot_threshold in  deg
if ~exist('interp_factor','var')|| isempty(interp_factor)
    interp_factor=[2,3.333];
end
tran_threshold=tran_threshold*1000;
rot_threshold=rot_threshold*1000;

prefix_fname=strtok2(filename(fname),'.');
pro=fullfile(strtok2(fname,'.'),[prefix_fname,'.pro']);

%[dsb,lin,par,sl,ushSet,timeStamp,freePara,dfn,lin_dfn,par_dfn,set_dfn,dref_dfn,d_ps] = readMeasDat(fsb,inf,0,true);
[dataStruct] = readMeasDat_savemat(fname,Inf);

Set=dataStruct.Set;
Line=dataStruct.Line;
Data=dataStruct.Data;
Rep=dataStruct.Repetition;

nch=32;
Set=Set(1:nch:end);
Line=Line(1:nch:end);
Rep=Rep(1:nch:end);
icePara=dataStruct.icePara(1:nch:end,:);

ncol=size(Data,1);
Data=reshape(Data,[ncol,nch,size(Data,2)/nch]);

nlin=double(max(Line))+1;
nrep=double(max(Rep))+1;
nset=double(max(Set))+1;

%pat=readsPar(pro,'lAccelFactPE');

nreacq=size(Data,3)-nlin*nset*nrep;
dro=reshape(Data(:,:,1:nlin*nset*nrep),[ncol,nch,nset,nlin,nrep]);
dreacq=reshape(Data(:,:,nlin*nset*nrep+1:end),[ncol,nch,nset,nreacq/nset]);
dro_rep=dro*0; %repeated acquisition

clear Data dataStruct;

PosCounter=1;  % 1 based; head positions
m_tran=0;
m_rot=0;
cReacquiredLine=0;

m_tran_orig=NaN(nlin,nrep); %the tran motion in the original (i.e. not reacquired) k-space
m_rot_orig=NaN(nlin,nrep);
m_tran_reacq=NaN(nlin,nrep); %the tran motion of the reacquired data (the least value if acquired multiple times)
m_rot_reacq=NaN(nlin,nrep); %the tran motion of the reacquired data

for i=1:length(Rep)
    
    tran=icePara(i,1);  %0.001 mm
    rot=icePara(i,2);   % 0.001 deg
    % motion parameter is changed at the end after kerenel run/ k-space data acquired
    if (tran~=m_tran(PosCounter) || rot ~=m_rot(PosCounter))
        PosCounter=PosCounter+1;
        m_rot(PosCounter)=rot;
        m_tran(PosCounter)=tran;
    end
    if i<=nlin*nset*nrep
        p2pos(Line(i)+1,Rep(i)+1)=PosCounter; %index to minimum m_rot and m_tran for the given k-space position.
        
    else
        if Set(i)==nset-1
            cReacquiredLine=cReacquiredLine+1;
            pReacquiredDataInfo(cReacquiredLine).rep=Rep(i);
            pReacquiredDataInfo(cReacquiredLine).lin=Line(i);
            pReacquiredDataInfo(cReacquiredLine).ipos=PosCounter;
        end
    end
end

for i=1:nlin
    for j=1:nrep
        [tran0,rot0]=get_motion(p2pos(i,j),m_tran,m_rot);
        m_tran_orig(i,j)=tran0;
        m_rot_orig(i,j)=rot0;
    end
end

for lReacqIndex=1:cReacquiredLine
    
    ipos=pReacquiredDataInfo(lReacqIndex).ipos;
    cRep=pReacquiredDataInfo(lReacqIndex).rep;
    cLin=pReacquiredDataInfo(lReacqIndex).lin;
    
    [tran,rot]=get_motion(ipos,m_tran,m_rot);
    
    if isnan(m_tran_reacq(cLin+1,cRep+1))||...
            (m_tran_reacq(cLin+1,cRep+1)+m_rot_reacq(cLin+1,cRep+1))>tran+rot
        m_tran_reacq(cLin+1,cRep+1)=tran;
        m_rot_reacq(cLin+1,cRep+1)=rot;
        dro_rep(:,:,:,cLin+1,cRep+1)=dreacq(:,:,:,lReacqIndex);
    end
    
end

m_to_reacq=m_tran_orig>tran_threshold | m_rot_orig>rot_threshold; %only data within the mask will be replaced, unlike m_reacq=~isnan(m_tran_reacqd);
%m_less_thr_rplc= m_to_reacq & (m_tran_reacq<=tran_threshold | m_rot_reacq<=rot_threshold);  %maybe | instead of &?
m_less_thr_rplc= m_to_reacq & (m_tran_reacq<=tran_threshold & m_rot_reacq<=rot_threshold);% changed to & on 11/10/2021


m_reacq=~isnan(m_tran_reacq);

m_to_reacq2=0*m_to_reacq;
m_less_thr_rplc2=0*m_less_thr_rplc;
m_reacq2=0*m_reacq;

m_to_reacq2(:,irep)=m_to_reacq(:,irep);
m_less_thr_rplc2(:,irep)=m_less_thr_rplc(:,irep);
m_reacq2(:,irep)=m_reacq(:,irep);

perc_to_reacq=sum(m_to_reacq2(:))/nlin/length(irep)*100;
perc_reacq_noneed=sum(m_reacq2(:)&~m_to_reacq2(:))/nlin/length(irep)*100;
perc_reacq_need=sum(m_reacq2(:)&m_to_reacq2(:))/nlin/length(irep)*100;
perc_less_thr_rplc=sum(m_less_thr_rplc2(:))/nlin/length(irep)*100;
%m_rplc=m_to_reacq2 & ((m_tran_reacq+m_rot_reacq)<(m_tran_orig+m_rot_orig) | m_less_thr_rplc2)&m_reacq2;
m_rplc=m_to_reacq2 & ((m_tran_reacq+m_rot_reacq)<(m_tran_orig+m_rot_orig) & m_less_thr_rplc2)&m_reacq2; % changed to & on 11/10/2021
perc_rplc = sum(m_rplc(:))/nlin/length(irep)*100;

dro2=dro;

m_rplc=shiftdim(m_rplc,-3);
sz=size(dro);
sz(4:5)=1;
dro2(repmat(m_rplc,sz))=dro_rep(repmat(m_rplc,sz));

dro=permute(dro,[1,4,2,3,5]);
dro=reshape(dro,[ncol,nlin,1,nch,nset,nrep]);

if ~exist(['Mag_',prefix,'.nii.gz'],'file')
  kdata2image(dro(:,:,:,:,:,irep),interp_factor,pro,prefix); 
end


tran_thr=num2str(tran_threshold/1000); %7/22/2021
rot_thr=num2str(rot_threshold/1000);
tran_thr=strrep(tran_thr,'.','p');
rot_thr=strrep(rot_thr,'.','p');

dro2=permute(dro2,[1,4,2,3,5]);
dro2=reshape(dro2,[ncol,nlin,1,nch,nset,nrep]);
if sum(m_rplc(:))>0
  kdata2image(dro2(:,:,:,:,:,irep),interp_factor,pro,[prefix,'_MCthr_',tran_thr,'_',rot_thr]);
end


save([prefix,'_perc_thr_',tran_thr,'_',rot_thr],'perc_to_reacq','perc_reacq_noneed','perc_reacq_need','perc_less_thr_rplc','perc_rplc',...
    'm_rot_reacq','m_tran_reacq','m_rot_orig','m_tran_orig','m_tran','m_rot','m_rplc','p2pos','irep');


function kdata2image(kdata,interp_factor,pro,prefix)
% kdata should be nro*npe*nsl*nch*nvenc*nt; image domain for the first dim

kdata= fft1c(kdata,1);


im = do_recon(kdata);

orient=proDimCenter(pro,round(size(im(:,:,1,1)).*interp_factor));


im2=do_interp_dim12_local(im,interp_factor);


save_nii(mean(im2,6),[prefix,'_mean'],kdata(:,1:32),orient);
save_nii(im2,prefix,kdata(:,1:32),orient);



function [tran,rot]= get_motion(ipos,m_tran,m_rot)


tran=0;
rot=0;

for i=ipos:ipos+2
    
    
    if (i<1) || i>length(m_tran)
        continue;
    end
    
    if (tran<m_tran(i))
        tran=m_tran(i);
    end
    if (rot<m_rot(i))
        rot=m_rot(i);
    end
end



%%


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

%
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

