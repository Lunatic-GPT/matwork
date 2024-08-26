function recon_fl_fq_pmoco(fname,interp_factor,prefix,irep)
%irep: the number of repetitions used for average calculation.
%7/12/2021: change the indexing in get_motion
if ~exist('interp_factor','var')|| isempty(interp_factor)
    interp_factor=[2,3.333];
end

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

nlin=max(Line)+1;
nrep=max(Rep)+1;
nset=max(Set)+1;

if ~exist('irep','var')
 irep=1:nrep;
end

pat=readsPar(pro,'lAccelFactPE');

nreacq=size(Data,3)-nlin*nset*nrep;
dro=reshape(Data(:,:,1:nlin*nset*nrep),[ncol,nch,nset,nlin,nrep]);
dreacq=reshape(Data(:,:,nlin*nset*nrep+1:end),[ncol,nch,nset,nreacq/nset]);

clear Data dataStruct;
    
PosCounter=1;  % 1 based; head positions
m_tran=0;
m_rot=0;
cReacquiredLine=0;
for i=1:length(Rep)
    
   	tran=icePara(i,1);
        rot=icePara(i,2);
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

fprintf('Total PosCounter = %d; total PosCounter without reacq = %d\n',PosCounter,max(p2pos(:)));
fprintf('Total k-space lines (no reacq) = %d; reacquired lines = %d\n',nlin*nrep,cReacquiredLine);

fprintf('tran/rot\n');
for i=1:length(m_tran)
  fprintf('%d/%d; ',m_tran(i),m_rot(i));
  if mod(i,10)==0
    fprintf('\n');
  end
end

dro2=dro;

replaceInfo=[];
replaceN=0;
for lReacqIndex=1:cReacquiredLine
		
    ipos=pReacquiredDataInfo(lReacqIndex).ipos;
    cRep=pReacquiredDataInfo(lReacqIndex).rep;
    cLin=pReacquiredDataInfo(lReacqIndex).lin;
	
    [tran,rot]=get_motion(ipos,m_tran,m_rot);
    
    [tran0,rot0]=get_motion(p2pos(cLin+1,cRep+1),m_tran,m_rot);
	
    if (tran0+rot0>tran+rot)	
         replaceN=replaceN+1;
         replaceInfo(replaceN,:)=[lReacqIndex,tran,rot,cLin,cRep,tran0,rot0];
          	
        dro2(:,:,:,cLin+1,cRep+1)=dreacq(:,:,:,lReacqIndex);	
	    p2pos(cLin+1,cRep+1)=ipos;
    end
	
end

%fprintf('%d/%d replaced/total acquired lines:\n lReacqIndex,tran,rot,cLin,cRep,tran0,rot0\n',...
%     size(replaceInfo,1),cReacquiredLine);
%disp(replaceInfo);

dro=permute(dro,[1,4,2,3,5]);
dro=reshape(dro,[ncol,nlin,1,nch,nset,nrep]);
if ~exist(['Mag_',prefix,'.nii.gz'],'file')
 kdata2image(dro(:,:,:,:,:,irep),interp_factor,pro,prefix);
end

if replaceN>0
 dro2=permute(dro2,[1,4,2,3,5]);
 dro2=reshape(dro2,[ncol,nlin,1,nch,nset,nrep]);
 kdata2image(dro2(:,:,:,:,:,irep),interp_factor,pro,[prefix,'_MC']);
end
save([prefix,'_motionPar'],'p2pos','m_tran','m_rot','cReacquiredLine','replaceN');


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

for i=ipos:ipos+2 %was ipos-1:ipos+1 before 7/12/2021
    
    
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

