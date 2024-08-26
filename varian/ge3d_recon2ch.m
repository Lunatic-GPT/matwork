function img =ge3d_recon2ch(fid_prefix)
% tbgemsShaper(fid_dir,format)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% format: output format: 'a': afni format. 's': sdt format.

fid_dir=[fid_prefix,'.fid'];

prefix = fid_prefix;
tic;
z = read_fid(fullfile(fid_dir,'fid'));
disp(toc);
arrvar=readPar(fid_dir,'array');
if length(arrvar)==2
    ne=1;
else
    arrval=readPar(fid_dir,arrvar(2:end-1));
    ne=length(arrval);
end

z=reshape(z,[size(z,1),size(z,2),2*ne,size(z,3)/2/ne]);
z=permute(z,[1,2,4,3]);
disp(toc);
z = squeeze(z);
%z = dcCorr(z);
disp(toc);
z = circshift(z,[size(z,1)/2,size(z,2)/2,size(z,3)/2,0]);
disp(toc);
z=fft(z,[],1);
disp(toc);
z=fft(z,[],2);
disp(toc);
z=fft(z,[],3);
disp(toc);

z = circshift(z,[size(z,1)/2,size(z,2)/2,size(z,3)/2,0]);
disp(toc);
lro = readPar(fid_dir,'lro');%in cm
lpe = readPar(fid_dir,'lpe');%in cm
lpe2 = readPar(fid_dir,'lpe2');%in mm

disp(toc);
    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    ppe2=readPar(fid_dir,'ppe2'); %in cm
    sz=size(z);
    delta = [lro*10/sz(1),lpe*10/sz(2),lpe2*10/sz(3)];
        
    orig_delta(1,:) = [pro,ppe,ppe2]*10-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    orig_delta(2,:) = delta;
    
    
 %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);

    if strcmp(orient,'sag');
       z=permute(z,[3,2,1,4]);
       orig_delta=orig_delta(:,[3,2,1]);
    end

    
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    
    aimg = sqrt(abs(z(:,:,:,1:2:end-1)).^2+abs(z(:,:,:,2:2:end)).^2);
    mag=cat(5,abs(z(:,:,:,1:2:end-1)),abs(z(:,:,:,2:2:end)),aimg);
    clear aimg;
    mag=permute(mag,[1,2,3,5,4]);
    szm=size(mag);
    mag=reshape(mag,[szm(1:3),szm(4)*szm(5)]);
    
    disp(toc);
    write_afni(mag,info,[prefix,'_mag']);
    disp(toc);
    clear mag;
    write_afni(angle(z),info,[prefix,'_ph']);
    
disp(toc);
    
    
%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');
if nargout ==0
 fprintf('ge3d_recon finished. Image dimensions:%d*%d*%d\n',size(z,1),size(z,2),size(z,3));
 img = '';
end

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
% ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave,[size(z_tmp,1),1,1,1]);
        

        
        
