function outname=do_interp_dim12(d,interp,prefix,nois,orient,do_scale,roi)
%do_interp_dim12(d,interp,mid,nois,suffix,do_scale,roi)
% im should be nro*npe*nsl*nch*nvenc*nt

% nois: nro*nch; can be empty
% phase_Only: only calculate the phase images, no mag or PC.
% 9/27/2020:
% now save in nii format
% deleted phase_only parameter;
% delete mid; now use prefix;
% add orient parameter;


if ~exist('do_scale','var')
    do_scale=true;
end

if ~exist('nois','var')
    nois='';
end

if ~exist('suffix','var')
    suffix='';
end

if ~exist('roi','var')
    roi=[];
end

d=ri(d);

nt=size(d,6);
max_mag=zeros(1,nt);

if max(roi(:))==1
    roi=clusterize2(roi);
end

sz=size(d);

sz(1:2)=round(sz(1:2).*interp);

if do_scale
    mag=zeros(sz(1),sz(2),sz(3),sz(5),nt,'uint16');
    pc=zeros(sz(1),sz(2),sz(3),sz(5)-1,nt,'int16');
    ph=zeros(sz(1),sz(2),sz(3),sz(5),nt,'int16');
else
    
    mag=zeros(sz(1),sz(2),sz(3),sz(5),nt,'single');
    pc=zeros(sz(1),sz(2),sz(3),sz(5)-1,nt,'single');
    ph=zeros(sz(1),sz(2),sz(3),sz(5),nt,'single');
    
    mag_4ph=zeros(sz(1),sz(2),sz(3),sz(5),nt,'single');
    
end

for l=1:nt
    d2=zeros(sz(1:5),'single');
    if any(interp>1)
        for j=1:sz(4)
            for k=1:sz(5)
                
                fd=single(fft2c(d(:,:,:,j,k,l)));
                tmp=zpad(fd,sz(1),sz(2),sz(3));
                
                d2(:,:,:,j,k)=ifft2c(tmp);
            end
        end
    else
        d2=d;
    end
    
    %     [mag(:,:,:,:,l),pc(:,:,:,:,l),max_mag(l)]=coil_combine_PC(d2,nois,false,do_scale);
    %     [tmp,ph(:,:,:,:,l)]=coil_combine(d2,nois,false);
    %

        [mag(:,:,:,:,l),pc(:,:,:,:,l),max_mag(l)]=coil_combine_PC(d2,nois,false,do_scale);
    
    for iroi=1:max(roi(:))
        [tmp,ph(:,:,:,:,l+(iroi-1)*nt)]=coil_combine(d2,nois,false,true,do_scale,roi==iroi);
    end
    
end

if do_scale
    for l=1:nt
        tmp=mag(:,:,:,:,l);
        mag(:,:,:,:,l)=uint16(single(tmp)*max_mag(l)/max(max_mag));
    end
end

%calculate pc_mean and mag_mean
%{
if nt>1
    mag_mean=mean(mag,5);
    dmean=zeros(sz(1),sz(2),sz(3),sz(5)-1,'single');
    for i=1:nt
        dmean=dmean+repmat(single(mag(:,:,:,1,i)),[1,1,1,sz(5)-1]).*single(mag(:,:,:,2:end,i)).*exp(single(pc(:,:,:,:,i))/18000*1i*pi);
    end
    pc_mean=int16(angle(dmean)*18000/pi);
end
%}
nro=size(mag,1);
nsl=size(mag,3);
nvenc=size(mag,4);

%mkdir(prefix);
%cur_dir=cd(prefix);
    
    mag=reshape(mag,[nro,size(mag,2),nsl,nvenc*nt]);
    pc=reshape(pc,[nro,size(pc,2),nsl,nt]);
    save_nii_local(mag,orient,['Mag_',prefix]);
    save_nii_local(pc,orient,['Phase_',prefix]);
        
%cd(cur_dir);

function save_nii_local(d,orient,prefix)

nii=make_nii(d);
nii=nii_from_orient(nii,orient);
   
save_untouch_niigz(nii,prefix);



%{
prefix=['recon_',mid,'_'];

mkdir([prefix,'interp',interps(2:end)]);
cur_dir=cd([prefix,'interp',interps(2:end)]);

if ~phase_Only
if nvenc==2
    
    mag=reshape(mag,[nro,size(mag,2),nsl,nvenc*nt]);
    pc=reshape(pc,[nro,size(pc,2),nsl,nt]);
    
    if nt>1
        d=mag;
      save(sprintf('Mag_%s_nt%d.mat',mid,nt),'d');
      d=pc;
      save(sprintf('Phase_%s_nt%d.mat',mid,nt),'d');
      outname=dir(sprintf('Mag_%s_nt%d.mat',mid,nt));
      outname(end+1)=dir(sprintf('Phase_%s_nt%d.mat',mid,nt));
    else
        d=mag;
        save(sprintf('Mag_%s.mat',mid),'d');
        d=pc;
        save(sprintf('Phase_%s.mat',mid),'d');
         outname=dir(sprintf('Mag_%s.mat',mid));
      outname(end+1)=dir(sprintf('Phase_%s.mat',mid));
    end
    
else
    mag=reshape(mag,[nro,size(mag,2),nsl,nvenc*nt]);
    if nt==1
      d=mag;
      save(sprintf('Mag_%s.mat',mid),'d');
      outname=dir(sprintf('Mag_%s.mat',mid));
    else
       d=mag;
       save(sprintf('Mag_%s_nt%d.mat',mid,nt),'d');
       outname=dir(sprintf('Mag_%s_nt%d.mat',mid,nt));
       
    end
  %  save(sprintf('Mag_4ph_%s.mat',mid),'mag_4ph');
    for i=1:nvenc-1
        pctmp=reshape(pc(:,:,:,i,:),[nro,size(pc,2),nsl,nt]);
        if nt==1
          d=pctmp;
          save(sprintf('Phase_%s_VENC%d',mid,i),'d');
          outname(end+1)=dir(sprintf('Phase_%s_VENC%d',mid,i));
        else
          d=pctmp;
          save(sprintf('Phase_%s_VENC%d_nt%d',mid,i,nt),'d');
          outname(end+1)=dir(sprintf('Phase_%s_VENC%d_nt%d',mid,i,nt));
        end
        
    end
    
end
end
if ~isempty(roi)
    d=ph;
    save(sprintf('PhaseImage_%s',mid),'d');
    outname(end+1)=dir(sprintf('PhaseImage_%s',mid));
end

if ~phase_Only
    if nt>1
        d=pc_mean;
        save(sprintf('Phase_%s_mean',mid),'d');
        d=mag_mean;
        save(sprintf('Mag_%s_mean.mat',mid),'d');
        outname(end+1)=dir(sprintf('Phase_%s_mean.mat',mid));
        outname(end+1)=dir(sprintf('Mag_%s_mean.mat',mid));
    end
end


%}
