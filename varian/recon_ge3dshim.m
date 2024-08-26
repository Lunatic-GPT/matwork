function [img,ph] =recon_ge3dshim(fid_dir,format)
% recon_ge3dshim(fid_dir,format)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% format: output format: 'a': afni format. 's': sdt format.

if ~exist('format','var')
    format = 'a';
end

prefix = fid_dir(1:end-4);

z = read_fid(fullfile(fid_dir,'fid'));
z = squeeze(z);
%z = dcCorr(z);

%z= circshift(z,[-6,0]);  %% phase correction.
z_shft = fftshift(z);
fz = fftn(z_shft);
fz2 = fftshift(fz);

img = abs(fz2);
ph = angle(fz2);
ph(isnan(ph))=0;
if isempty(prefix)
    prefix = 'a';
end

lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
lpe2 = readPar(fid_dir,'lpe2'); %in cm
  
img=flipdim(img,3);
ph=flipdim(ph,3);

img=flipdim(img,2);
ph=flipdim(ph,2);


orient=readPar(fid_dir,'orient');
    
sz=size(img);

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    ppe2=readPar(fid_dir,'ppe2'); %in cm
    
    delta = [lro*10/sz(1),lpe*10/sz(2),lpe2*10/sz(3)];
    
    orig_delta(1,:) = [pro,ppe,ppe2]*10-delta.*(sz/2-[0.5,0.5,0.5]);
    orig_delta(2,:) = delta;
    
    [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
    ph = reorient_data(ph,[],orient(2:end-1));
    
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    info.HISTORY_NOTE = 'recon_ge3dshim';
if format == 'a'
    write_afni(img,info,prefix);
 %   write_afni(ph,info,[prefix,'_ph']);
else 
   
    ad = abs(orig_delta(2,:)).*[size(img,1),size(img,2),1];
     writeSdt2(img,prefix,ad(1),ad(2),ad(3),1);
end

%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');
if nargout ==0
 fprintf('recon_ge3dshim shaper finished. Image dimensions:%d*%d*%d\n',size(img,1),size(img,2),size(img,3));
 img = '';
end

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
 ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave3,[size(z_tmp,1),size(z_tmp,2),1,1]);
        

