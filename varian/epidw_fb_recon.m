function [img,ph] =epidw_fb_recon(fid_dir,dopartialFT,keep_refscan,nophase)
% epiShaper_1dir_halfFT(fid_dir,halfFT,keep_refscan,nophase)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% keep_refscan: keep the image given by ref_pos+1.  Default: false.
% nophase: do not save phase image.default true.
if ~strcmp(computer,'PCWIN')
    if exist('/data/users/xiaopeng/vnmrsys','dir')
      root='/data/users/xiaopeng/vnmrsys';
    elseif exist('/home/xiaopeng/vesta/vnmrsys','dir');
      root='/home/xiaopeng/vesta/vnmrsys';
    else
        root='/home/xiaopeng/labhome/vnmrsys';
    end
     
else
    root = 'z:/home/xiaopeng/vnmrsys';
end

if ~exist('nophase','var')
    nophase = true;
end

if ~exist('dopartialFT','var')
    dopartialFT = true;
end

if ~exist('keep_refscan','var')
    keep_refscan = true;
end

if ~exist(fid_dir,'dir')
  fid_dir=[fid_dir,'.fid'];
end

[pathstr,prefix]=fileparts(fid_dir);

if exist(fullfile(fid_dir,'fid.orig'),'file')
    fid_file = 'fid.orig';
    procpar = 'procpar.orig';
else
    fid_file = 'fid';
    procpar = 'procpar';
end

z_tmp = read_fid(fullfile(fid_dir,fid_file));
ndim=ndims(z_tmp);

for i=1:size(z_tmp,3)
    z_tmp2=z_tmp(:,:,i);
    if ~any(z_tmp2(:)>0)
        z_tmp=z_tmp(:,:,1:i-1);
        break;
    end
end
        

ns =  readPar(fid_dir,'ns');
nseg = readPar(fid_dir,'nseg');
nv =  readPar(fid_dir,'nv');
npe=size(z_tmp,2)/ns/nseg;

z_tmp = reshape(z_tmp,[size(z_tmp,1),npe,ns*nseg,size(z_tmp,3)]);


drank=[3,0];
if exist([prefix,'+orig.HEAD'],'file')
    [err,info] = BrikInfo([prefix,'+orig.HEAD']);
    drank = info.DATASET_RANK;
    if drank(2)>size(z_tmp,4)
        fprintf('Delete old file: %s+orig\n',prefix);
        delete([prefix,'+orig*']);
        drank(2)=0;
    elseif drank(2)==size(z_tmp,4)
  %      fprintf('No new data acquired\n');
         img=[];
        return;
    else    
     
    % fprintf('Add %d briks to %s+orig\n',size(z_tmp,4),prefix);
    end
end

%z = convertTraces(z_tmp);  % reverse the direction of even number traces.
z = dcCorr(z_tmp);
image = readPar(fid_dir,'image');
refpos=find(image==0);
if refpos>size(z,4)
    error('reference scan not yet acquired');
end

z = fidPhaseCorr(z,refpos,keep_refscan);

z = z(:,:,:,drank(2)+1:end);

%z= circshift(z,[-6,0]);  %% phase correction.

% segmentation re-order
z=reshape(z,[size(z,1),npe,ns,nseg,size(z,4)]);
z=permute(z,[1,2,4,3,5]);
z=reshape(z,[size(z,1),nseg*npe,ns,size(z,5)]);
petable = readPar(fid_dir,'petable');
etl2=readPar(fid_dir,'etl2');
fract_ky=readPar(fid_dir,'fract_ky');
i_skip=size(z,2)-etl2-1;


z(:,i_skip:i_skip+1,:,:)=[];

z2=tabc(z,fullfile(root,'tablib',petable(2:end-1)),'t1');

sz=size(z2);
%sz(2)=nv/2;
sz(2)=size(z2,2)/2;
sz2=sz;
sz2(4)=sz2(4)*2;
if dopartialFT
    sz2(2)=nv/2;
end    
img2 = zeros(sz2);
ph2=zeros(sz2);

for iimg=1:2
z = z2(:,iimg:2:end,:,:);  % assuming linear order.
if ~dopartialFT
  [img,ph]=normalFT(z);
else
  [img,ph]=partialFT(z,fract_ky/2);
end
lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
thk = readPar(fid_dir,'thk');%in mm
total_tr = readPar(fid_dir,'tr');
pss=readPar(fid_dir,'pss'); %in cm
orient=readPar(fid_dir,'orient');
img = slice_reorder(img,pss);    
ph = slice_reorder(ph,pss);    
img = flipdim(img,1);  %make it left to right
ph = flipdim(ph,1);   

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    sz=size(img);
    delta = [lro*10/sz(1),lpe*10/sz(2),thk];  % thickness in mm
    if length(pss)>1
        pss_sort = sort(pss);
        delta(3) = (pss_sort(2)-pss_sort(1))*10;  %pss in cm 
    end
    
    orig_delta(1,:) = [pro,ppe,min(pss)]*10-delta.*([sz(1:2)/2,0]-[0.5,0.5,0]);
    orig_delta(2,:) = delta;
    [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
    ph = reorient_data(ph,[],orient(2:end-1));
    img2(:,:,:,(iimg-1)*sz(4)+(1:sz(4)))=img;
    ph2(:,:,:,(iimg-1)*sz(4)+(1:sz(4)))=ph;
end

    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);   % change from varian convention to afni convention.
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    info.TAXIS_FLOATS = [0,total_tr/1000,0,0,0];

    if exist([prefix,'+orig.HEAD'],'file');
       write_afni(img,info,'epiShaper_tmp');
        cmd = sprintf('3dTcat -glueto %s+orig %s+orig',prefix,'epiShaper_tmp');
      unix(cmd);
      delete epiShaper_tmp*
    else
        write_afni(img2,info,prefix);
    end
      
    if ~nophase
     write_afni(ph2,info,[prefix,'_ph']);
    end


if nargout ==0
 fprintf('epiShaper finished. Image dimensions:%d*%d*%d*%d\n',size(img));
 img = [];
end

%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');


function z = convertTraces(z_tmp)
 % reverse the direction of even number traces
 z = zeros(size(z_tmp));
 
 z(:,1:2:end-1,:,:) = z_tmp(:,1:2:end-1,:,:);
 z(1:end,2:2:end,:,:) = z_tmp(end:-1:1,2:2:end,:,:);

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
 ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave3,[size(z_tmp,1),size(z_tmp,2),1,1]);
        
 
function z=fidPhaseCorr(z,refpos,keep_refscan)
  %
 NUM_DC_SUM = 8;
 
        
        zref = z(:,:,:,refpos);
        sz = size(z);
        dc_pts = [1:NUM_DC_SUM/2,sz(1)-NUM_DC_SUM/2+1:sz(1)];
        zref = zref-repmat(mean(zref(dc_pts,:,:),1),[sz(1),1,1]);
        
        z=fftshift(z,1);  %shift 
        zref=fftshift(zref,1); %shift
        
        fzref = fft(zref,[],1);
        fzref = fzref./abs(fzref);  % normalize
        
        
        
        sz_rep = sz;
        sz_rep(1:end-1)=1;
        fzref = repmat(fzref,sz_rep);
        
        fz = fft(z,[],1);
        fz = fz.*conj(fzref);
        
        z=ifft(fz,[],1);
        
        z=fftshift(z,1);  %shift back
        if ~keep_refscan
          z(:,:,:,refpos)=[];  %remove the reference scan.
        end
        
function [img,ph]=normalFT(z)

 sz=size(z);
 
        img=zeros(sz);
        ph=zeros(sz);
  for i=1:size(z,3)
    for j=1:size(z,4)
        z_shft = fftshift(z(:,:,i,j));
        fz = fft2(z_shft);
        fz2 = fftshift(fz);
       img(:,:,i,j) = abs(fz2);
       ph(:,:,i,j) = angle(fz2);
    end
  end
ph(isnan(ph))=0;
