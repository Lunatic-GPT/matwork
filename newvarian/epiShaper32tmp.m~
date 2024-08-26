function epiShaper32tm(fid_dir,triple_ref,nophase)
% epiShaper32(fid_dir,nophase)
% nophase: do not save phase image.default true.


if ~exist('nophase','var')
    nophase = true;
end

if ~exist(fid_dir,'dir')
  fid_dir=[fid_dir,'.fid'];
end

if ~exist('triple_ref','var')
    triple_ref=true;
end

[pathstr,prefix]=fileparts(fid_dir);

z_tmp = read_fid(fid_dir);
nnav=readPar(fid_dir,'nnav');
nread=readPar(fid_dir,'nread');
nphase=readPar(fid_dir,'nphase');
nseg=readPar(fid_dir,'nseg');
kzero=readPar(fid_dir,'kzero');
ns = readPar(fid_dir,'ns');
im = readPar(fid_dir,'image');
os = readPar(fid_dir,'oversample');  %oversample when lro is not equal to 0.
nav_type=readPar(fid_dir,'nav_type');
na=length(im);
z=reshape(z_tmp,[nread*os/2,nphase/nseg+nnav+1,ns,nseg,na]);

pro=readPar(fid_dir,'pro');
ppe=readPar(fid_dir,'ppe');
lro=readPar(fid_dir,'lro');
lpe=readPar(fid_dir,'lpe');


z = convertTraces(z);  % reverse the direction of even number traces.
z=fov_shift_kspace(z,[-pro,-ppe],[lro,lpe]);
z=z(:,[1:nnav,nnav+2:end],:,:,:);
z=permute(z,[1,2,4,3,5]);


z=circshift(z,[size(z,1)/2,0,0,0,0]);

fxz=fft(z,[],1);

fxz=circshift(fxz,[size(z,1)/2,0,0,0,0]);


%z=reshape(z,[size(z,1),nphase+nseg*nnav,ns,na]);


%z = dcCorr(z);
epiref_type=readPar(fid_dir,'epiref_type');

if strcmp(epiref_type,'pointwise') || ~triple_ref
 z2 = fidPhaseCorr(z(:,:,:,im>0),z(:,:,:,im==0));

  % z2=navcorr(z2,nnav,nav_type);
   
   z2=permute(z2,[1,3,2,4,5]);
   
   z2=reshape(z2,[size(z2,1),size(z2,2)*size(z2,3),size(z2,4),size(z2,5)]);
   
   z2=circshift(z2,[0,nphase/2+1,0,0]);
   
   z2=fft(z2,[],2);
   z2=circshift(z2,[0,nphase/2,0,0]);
   
else
    
   zr1 =fxz(:,:,:,:,im==0);
   zr2 = fxz(:,:,:,:,im==-2);
   zrp = fxz(:,:,:,:,im==-1);
   zdata=fxz(:,:,:,:,im==1);
   if size(zrp,5)==1
       zrp=repmat(zrp,[1,1,1,1,size(zdata,5)]);
   end
   
   
   z2=fidPhaseCorr(zdata,zr1);
   zrp2=fidPhaseCorr(zrp,zr2);
   
   %nav_type='"off"';
   z2=navcorr(z2,nnav,nav_type);
   zrp2=navcorr(zrp2,nnav,nav_type);
   
   z2=permute(z2,[1,3,2,4,5]);
   zrp2=permute(zrp2,[1,3,2,4,5]);
   
   z2=reshape(z2,[size(z2,1),size(z2,2)*size(z2,3),size(z2,4),size(z2,5)]);
   zrp2=reshape(zrp2,[size(zrp2,1),size(zrp2,2)*size(zrp2,3),size(zrp2,4),size(zrp2,5)]);
   

   zrp2=circshift(zrp2,[0,nphase/2+1,0,0]);
   zrp2=fft(zrp2,[],2);
   zrp2=circshift(zrp2,[0,nphase/2,0,0]);
   zrp2(2:end,:,:,:)=zrp2(end:-1:2,:,:,:);
   zrp2(1,:,:,:)=0;
   
   z2=circshift(z2,[0,nphase/2+1,0,0]);
   z2=fft(z2,[],2);
   z2=circshift(z2,[0,nphase/2,0,0]);
   

   scl=z2./repmat(z2(:,:,:,1),[1,1,1,size(z2,4)]);
   %scl=abs(scl);
   tmp=z2(:,:,:,1);
   mask=abs(tmp)>abs(max(tmp(:)))*0.01;
   mask=repmat(mask,[1,1,1,size(z2,4)]);
   z2(mask)=z2(mask)+zrp2(mask).*scl(mask);

end

% segmentation re-order

pescheme = readPar(fid_dir,'pescheme');

if ~isempty(pescheme) && ~strcmp(pescheme(2),'l')
  error('Only linear scheme has been implemented');
end


pss=readPar(fid_dir,'pss'); %in cm

z2 = slice_reorder(z2,pss); 

orient=readPar(fid_dir,'orient');
orient=orient(2:end-1);


switch orient
    case 'trans90'
     z2=flipdim(z2,2);  %up down
     z2=flipdim(z2,1); % left right  %left on the left side in stimulate.  
     afni=z2;
    case 'trans'
     z2 = permute(z2,[2,1,3,4]);   %do I really need this? yes
     afni=z2;
    case 'sag90'
      z2=flipdim(z2,1);
      z2=flipdim(z2,2);  
      afni=permute(b,[3,1,2,4]);
    case 'sag'
      z2=permute(z2,[2,1,3,4]);   
      z2=flipdim(z2,2);      %first dimesion is l-r in stimulate
      afni=permute(b,[3,1,2,4]);  
    otherwise
        
    warning('The directions may be wrong');
end


%img = flipdim(img,1);  %make it left to right
%ph = flipdim(ph,1);   
    
  writesdt4(abs(z2),[prefix,'_mag']);
  if ~nophase
    writesdt4(angle(z2),[prefix,'_ph']);
  end

 
thk = readPar(prefix,'thk');%in mm
tr = readPar(prefix,'tr');
%img = slice_reorder(img,pss);    

    pro=readPar(prefix,'pro'); %in cm
    ppe=readPar(prefix,'ppe'); %in cm
    sz=size(afni);
    delta = [lro*10/sz(1),lpe*10/sz(2),thk];  % thickness in mm
    if length(pss)>1
        pss_sort = sort(pss);
        delta(3) = (pss_sort(2)-pss_sort(1))*10;  %pss in cm 
    end
    
    
    delta([1,3])=-delta([1,3]);  %afni convention R2L, A2P, I2S
             
    switch orient
        
    %varian convention:
    %positive pro: axial P; axial90 R; coronal I; coronal90 R; sag I; sag90 P;
    %positive ppe: axial L; axial90 P; coronal L; coronal90 I; sag A; sag90 I;
    %positive pss: axial I; axial90 I; coronal A; coronal90 A; sag R; sag90 R;
    %     
        case 'trans'
            center = [-ppe,pro,min(pss)]*10; 
            center([1,3])=-center([1,3]);    
            orig = center-delta.*([sz(1:2)/2,0]-[0.5,0.5,0]);
            %[pro,ppe,min(pss)]*10;
        case 'trans90'
            center = [pro,ppe,min(pss)]*10;
            center([1,3])=-center([1,3]);    
            orig = center-delta.*([sz(1:2)/2,0]-[0.5,0.5,0]);
        case 'sag'
            center = [min(pss),-ppe,pro]*10;
            center([1,3])=-center([1,3]);    
            orig = center-delta.*([0,sz(2:3)/2]-[0,0.5,0.5]);
        case 'sag90'
            center = [min(pss),pro,ppe]*10;
            center([1,3])=-center([1,3]);    
            orig = center-delta.*([0,sz(2:3)/2]-[0,0.5,0.5]);
        otherwise
            
    end
          
    %    
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
     
    info.ORIGIN = orig;
    info.DELTA = delta;
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I
    write_afni(afni,info,d);
    
    
function z2=navcorr(z,nnav,nav_type)
 
    if strcmp(nav_type(2:end-1),'off')    
      z2=z(:,nnav+1:end,:,:,:);
    else
      ref=z(:,nnav,:,:,:);
      ref=repmat(ref,[1,size(z,2)-nnav,1,1,1]);
      
      corr=conj(ref)./abs(ref);
      corr(isnan(corr))=1;
      corr(isinf(corr))=1;
      z2=z(:,nnav+1:end,:,:,:).*corr;
    end
    
   
%seqcon = read_par(fullfile(fid_dir,'procpar'),'seqcon');
function img=normalFT(z)
  sz=size(z);
        img=zeros(sz);
 
  for i=1:size(z,3)
    for j=1:size(z,4)
       z_shft = fftshift(z(:,:,i,j),2);
     %z_shft = z(:,:,i,j);
     %z_shft=fftshift(z_shft);  
     fz = fft2(z_shft);
        
       fz2 = fftshift(fz);
       img(:,:,i,j) = fz2;
    end
  end
%ph(isnan(ph))=0;

function z = convertTraces(z_tmp)
 % reverse the direction of even number traces
 z = zeros(size(z_tmp));
 
 z(:,1:2:end-1,:,:,:) = z_tmp(:,1:2:end-1,:,:,:);
 z(1:end,2:2:end,:,:,:) = z_tmp(end:-1:1,2:2:end,:,:,:);

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
 ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave3,[size(z_tmp,1),size(z_tmp,2),1,1]);
        
 
function z=fidPhaseCorr(z,fzref)
  %
        sz = size(z);
        
      %  fzref = ifft(zref,[],1);
        fzref = fzref./abs(fzref);  % normalize
        sz_rep = sz;
        if ndims(z)>ndims(fzref)  
          sz_rep(1:end-1)=1;
          fzref = repmat(fzref,sz_rep);
        end
        
        z = z.*conj(fzref);
        
        
        
       
        
   function z=fidPhaseCorr2(z,zref)
  %
        sz = size(z);
        fz = fft(z,[],1);
        fzref = fft(zref,[],1);
        sz_rep = sz;
        if length(sz)==4
          sz_rep(1:end-1)=1;
          fzref = repmat(fzref,sz_rep);
        end
        H=mean(abs(fz(:)));
        ph=conj(fz).*fzref./abs(fz)./abs(fzref);
        ph(:,1:2:end-1,:,:)=1;
        %ph(abs(fz)<0.1*H & abs(fzref)<0.1*H)=1;
        
        
        fz = fz.*ph;
        
        z=ifft(fz,[],1);
        
        
        
        