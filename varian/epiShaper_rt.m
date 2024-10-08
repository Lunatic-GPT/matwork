function [img,n_proc] =epiShaper_rt(fid_dir,n_proc,keep_refscan,nophase)
% tbgemsShaper(fid_dir,keep_refscan,nophase)
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

if ~exist('keep_refscan','var')
    keep_refscan = true;
end

[pathstr,fid_name]=fileparts(fid_dir);

prefix = strtok(fid_name,'.');

if ~exist(fid_dir,'dir')
  fid_dir=[fid_dir,'.fid'];
end


fid_file = 'fid';

copyfile(fullfile(fid_dir,fid_file),'fid_temp');
z_tmp = read_fid('fid_temp');

for i=1:size(z_tmp,3)
    z_tmp2=z_tmp(:,:,i);
    if ~any(z_tmp2(:)>0)
        z_tmp=z_tmp(:,:,1:i-1);
        break;
    end
end

z_tmp=z_tmp(:,:,1:end-1); % do not process the last one.
        

ns =  readPar(fid_dir,'ns');
nseg = readPar(fid_dir,'nseg');
npe=size(z_tmp,2)/ns/nseg;

z_tmp = reshape(z_tmp,[size(z_tmp,1),npe,ns*nseg,size(z_tmp,3)]);


refpos = readPar(fid_dir,'ref_pos')+1;
if refpos>size(z_tmp,4)
    warning('reference scan not yet acquired');
    
    img=[];
    return;
end


z_ref=z_tmp(:,:,:,refpos);
z_tmp(:,:,:,refpos)=[];

    if n_proc>size(z_tmp,4)
        n_proc=0;
    elseif n_proc==size(z_tmp,4)
  %      fprintf('No new data acquired\n');
         img=[];
        return;
    else    
     
    % fprintf('Add %d briks to %s+orig\n',size(z_tmp,4),prefix);
    end

z=z_tmp(:,:,:,n_proc+1:end);
z=cat(4,z_ref,z);

z = convertTraces(z);  % reverse the direction of even number traces.
z = dcCorr(z);

z = fidPhaseCorr(z,1,false);

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

img = zeros(size(z));
ph=zeros(size(z));
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

lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
thk = readPar(fid_dir,'thk');%in mm
total_tr = readPar(fid_dir,'total_tr');
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
    n_proc=size(img,4)+n_proc;



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
        
        fzref = fft(zref,[],1);
        fzref = fzref./abs(fzref);  % normalize
        sz_rep = sz;
        sz_rep(1:end-1)=1;
        fzref = repmat(fzref,sz_rep);
        
        fz = fft(z,[],1);
        fz = fz.*conj(fzref);
        
        z=ifft(fz,[],1);
        if ~keep_refscan
          z(:,:,:,refpos)=[];  %remove the reference scan.
        end
        




