function [img,ph] =epiShaper(fid_dir,format,keep_refscan,nophase)
% tbgemsShaper(fid_dir,format,keep_refscan,nophase)
% format: 'a' or 's'
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% keep_refscan: keep the image given by ref_pos+1.  Default: false.
% nophase: do not save phase image.default true.

if ~exist('format','var')
    format='a';
end

if ~strcmp(computer,'PCWIN64')
    if exist('/data/users/xiaopeng/vnmrsys','dir')
      root='/data/users/xiaopeng/vnmrsys';
    elseif exist('/home/xiaopeng/vesta/vnmrsys','dir');
      root='/home/xiaopeng/vesta/vnmrsys';
    elseif strcmp(computer,'GLNXA64')
        root = '/home/xiaopeng/oldvnmrsys';
    else
        root='/home/xiaopeng/labhome/vnmrsys';
    end
     
else
    root = 'c:/labhome/vnmrsys';
end

if ~exist('nophase','var')
    nophase = true;
end

if ~exist('keep_refscan','var')
    keep_refscan = false;
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
    %    z_tmp=z_tmp(:,:,1:i-1);
        break;
    end
end
        

ns =  readPar(fid_dir,'ns');
nseg = readPar(fid_dir,'nseg');
npe=size(z_tmp,2)/ns/nseg;

z_tmp = reshape(z_tmp,[size(z_tmp,1),npe,ns*nseg,size(z_tmp,3)]);
drank=[3,0];
if 0
if exist([prefix,'+orig.HEAD'],'file')
    [err,info] = BrikInfo([prefix,'+orig.HEAD']);
    drank = info.DATASET_RANK;
    if drank(2)>size(z_tmp,4)
        fprintf('Delete old file: %s+orig\n',prefix);
        delete([prefix,'+orig*']);
        drank(2)=0;
    elseif drank(2)==size(z_tmp,4)  || (~keep_refscan && drank(2)==size(z_tmp,4)-1)
       fprintf('No new data acquired\n');
         img=[];
        return;
    else    
     
    % fprintf('Add %d briks to %s+orig\n',size(z_tmp,4),prefix);
    end
end
end
z = convertTraces(z_tmp);  % reverse the direction of even number traces.
z = dcCorr(z);
refpos = readPar(fid_dir,'ref_pos')+1;
if refpos>size(z,4)
    error('reference scan not yet acquired');
end

z = fidPhaseCorr(z,refpos,keep_refscan);

z = z(:,:,:,drank(2)+1:end);
shiftecho=readPar(fid_dir,'shiftecho');
negphases=readPar(fid_dir,'negphases');

if shiftecho=='y'
    z=circshift(z,[0,size(z,2)/2-negphases-1]);
end
%z= circshift(z,[-6,0]);  %% phase correction.

% segmentation re-order
z=reshape(z,[size(z,1),npe,ns,nseg,size(z,4)]);
z=permute(z,[1,2,4,3,5]);
z=reshape(z,[size(z,1),nseg*npe,ns,size(z,5)]);
petable = readPar(fid_dir,'petable');

z=tabc(z,fullfile(root,'tablib',petable(2:end-1)),'t1');
if ~shiftecho
[img,ph]=normalFT(z);
else
  [img,ph]=normalFT(z);
%[img,ph]=partialFT(z,negphases);
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
    
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);   % change from varian convention to afni convention.
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    info.TAXIS_FLOATS = [0,total_tr/1000,0,0,0];

  
    
    if exist([prefix,'+orig.HEAD'],'file');
       write_afni_sdt(img,info,'epiShaper_tmp',format);
        cmd = sprintf('3dTcat -glueto %s+orig %s+orig',prefix,'epiShaper_tmp');
      unix(cmd);
      delete epiShaper_tmp*
    else
        write_afni_sdt(img,info,prefix,format);
    end
      
    if ~nophase
     write_afni_sdt(ph,info,[prefix,'_ph'],format);
    end


if nargout ==0
 fprintf('epiShaper finished. Image dimensions:%d*%d*%d*%d\n',size(img));
 img = [];
end

%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');
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
        
   function z=fidPhaseCorr2(z,refpos,keep_refscan)
  %
        sz = size(z);
        
        z=fftshift(z,1);  %shift 
        fz = fft(z,[],1);
        
        fz=fftshift(fz,1);
        
        fzref = fz(:,:,:,refpos);
        
        afz=abs(fzref);
        ph=zeros(sz(1:3));
        for i=1:sz(3)
         ph(:,:,i)=phase2d(fzref(:,:,i),1);
        end
        ph_fit=zeros(sz(1:end-1));
        for i=1:size(fzref,2)
            for j=1:size(fzref,3)
                x=1:size(fzref,1);
                y=ph(:,i,j);
                ind=find(afz(:,i,j)>max(afz(:,i,j))*0.4);
                p=polyfit(ind,y(ind),1);
                ph_fit(:,i,j)=p(1)*x+p(2);
                
            end
        end
        sz_rep = sz;
        sz_rep(1:end-1)=1;
        ph_fit = repmat(ph_fit,sz_rep);
        fz = fz.*exp(-1i*ph_fit);
        
        fz=fftshift(fz,1);
        z=ifft(fz,[],1);
        z=fftshift(z,1);  %shift back
        if ~keep_refscan
          z(:,:,:,refpos)=[];  %remove the reference scan.
        end     

        
        
        