function img =ge3d_recon(fid_prefix,format)
% tbgemsShaper(fid_dir,format)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% format: output format: 'a': afni format. 's': sdt format.



if ~exist('format','var')
    format = 'a';
end

fid_dir=[fid_prefix,'.fid'];

prefix = fid_prefix;

z = read_fid(fullfile(fid_dir,'fid'));
z = squeeze(z);
z = dcCorr(z);

nv=readPar(fid_dir,'nv');
nv2=readPar(fid_dir,'nv2');
np=readPar(fid_dir,'np');

z=reshape(z,[np/2,nv,nv2]);

z_shft = fftshift(z);
fz = fftn(z_shft);
fz2 = fftshift(fz);

img = fz2;
if isempty(prefix)
    prefix = 'a';
end

lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
lpe2 = readPar(fid_dir,'lpe2');%in mm

img=flipdim(img,1);


sz=size(img);
if length(sz)==2
    sz(3)=1;
end

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    ppe2=readPar(fid_dir,'ppe2'); %in cm
    
    delta = [lro*10/sz(1),lpe*10/sz(2),lpe2*10/sz(3)];
        
    orig_delta(1,:) = [pro,ppe,ppe2]*10-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    orig_delta(2,:) = delta;
    
 %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);

if format == 'a'
    write_afni(abs(img),prefix,info);
    write_afni(angle(img),[prefix,'_ph'],info);
    
else 
    ad = abs(orig_delta(2,:)).*sz;
    writesdt4(abs(img),prefix);
end

%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');
if nargout ==0
 fprintf('ge3d_recon finished. Image dimensions:%d*%d*%d\n',size(img,1),size(img,2),size(img,3));
 img = '';
end

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
% ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave,[size(z_tmp,1),1,1,1]);
        

