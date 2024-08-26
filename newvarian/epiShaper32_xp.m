function epiShaper32_xp(fid_dir,triple_ref,format,nophase)
% epiShaper32(fid_dir,triple_ref,format,nophase)
% fid_dir: name of the fid directory.
% triple_ref: true for triple reference and false for single reference.
% nophase: do not save phase image.  default true.
% format: 'a' or 's', (afni or sdt). default 'a';
if ~exist('format','var')
    format='a';
end

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
nread=nread*os/2;
lro=lro*os;

z = convertTraces(z);  % reverse the direction of even number traces.
z(:,:,:,:,im==0|im==1)=fov_shift_kspace(z(:,:,:,:,im==0|im==1),[pro,0],[lro,lpe]);

z(:,:,:,:,im==-1|im==-2)=fov_shift_kspace(z(:,:,:,:,im==-1|im==-2),[-pro,0],[lro,lpe]);

z=z(:,[1:nnav,nnav+2:end],:,:,:);
z=permute(z,[1,2,4,3,5]);


z=circshift(z,[size(z,1)/2,0,0,0,0]);

fxz=fft(z,[],1);

fxz=circshift(fxz,[size(z,1)/2,0,0,0,0]);


%z=reshape(z,[size(z,1),nphase+nseg*nnav,ns,na]);


%z = dcCorr(z);
epiref_type=readPar(fid_dir,'epiref_type');
epiref_type=epiref_type(2:end-1);
if strcmp(epiref_type,'pointwise') || ~triple_ref
   z2 = fidPhaseCorr(fxz(:,:,:,:,im>0),fxz(:,:,:,:,im==0));
 %  z2=fxz;
   nav_type='"off"';
   z2=navcorr(z2,nnav,nav_type);
   
   z2=permute(z2,[1,3,2,4,5]);
   
   z2=reshape(z2,[size(z2,1),size(z2,2)*size(z2,3),size(z2,4),size(z2,5)]);
   
   z2=circshift(z2,[0,nphase/2+1,0,0]);
   if nseg==1
    z2=fov_shift_kspace(z2,[0,-ppe],[lro,lpe]);
   else
    z2=fov_shift_kspace(z2,[0,ppe],[lro,lpe]);
   end
   
   z2=fft(z2,[],2);
   z2=circshift(z2,[0,nphase/2,0,0]);

   
else
    
   zr1 =fxz(:,:,:,:,im==0);
   zr2 = fxz(:,:,:,:,im==-2);
   zrp = fxz(:,:,:,:,im==-1);
   zdata=fxz(:,:,:,:,im==1);
   if strcmp(epiref_type,'triple')
       zrp=repmat(zrp,[1,1,1,1,size(zdata,5)]);
   end
   
   
   z2=fidPhaseCorr(zdata,zr1);
   zrp2=fidPhaseCorr(zrp,zr2);
   
 %  nav_type='"off"';
   z2=navcorr(z2,nnav,nav_type);
   zrp2=navcorr(zrp2,nnav,nav_type);
   
   z2=permute(z2,[1,3,2,4,5]);
   zrp2=permute(zrp2,[1,3,2,4,5]);
   
   z2=reshape(z2,[size(z2,1),size(z2,2)*size(z2,3),size(z2,4),size(z2,5)]);
   zrp2=reshape(zrp2,[size(zrp2,1),size(zrp2,2)*size(zrp2,3),size(zrp2,4),size(zrp2,5)]);
   
   if nseg==1
     z2=fov_shift_kspace(z2,[0,-ppe],[lro,lpe]);
     zrp2=fov_shift_kspace(zrp2,[0,-ppe],[lro,lpe]);
   elseif nseg==4
     z2=fov_shift_kspace(z2,[0,ppe],[lro,lpe]);
     zrp2=fov_shift_kspace(zrp2,[0,ppe],[lro,lpe]);
   end
   
   zrp2=circshift(zrp2,[0,nphase/2+1,0,0]);
   zrp2=fft(zrp2,[],2);
   zrp2=circshift(zrp2,[0,nphase/2,0,0]);
   zrp2(2:end,:,:,:)=zrp2(end:-1:2,:,:,:);
   zrp2(1,:,:,:)=0;
   
   z2=circshift(z2,[0,nphase/2+1,0,0]);
   z2=fft(z2,[],2);
   z2=circshift(z2,[0,nphase/2,0,0]);
   
   if strcmp(epiref_type,'triple')
      scl=z2./repmat(z2(:,:,:,1),[1,1,1,size(z2,4)]);
      %scl=abs(scl);
      tmp=z2(:,:,:,1);
      mask=abs(tmp)>abs(max(tmp(:)))*0.01;
      mask=repmat(mask,[1,1,1,size(z2,4)]);
      %  mask=ones(size(z2));
      z2(mask>0)=0.5*(z2(mask>0)+zrp2(mask>0).*scl(mask>0));
   elseif strcmp(epiref_type,'fulltriple')
      z2=0.5*(z2+zrp2);
   end
end

pss=readPar(fid_dir,'pss'); %in cm

z2 = slice_reorder(z2,pss); 

 writesdt4(abs(z2),[prefix,'_mag']);
 if ~nophase
  writesdt4(angle(z2),[prefix,'_ph']);
 end

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
     %   fzref2=fzref(:,1,:,:,:);
      %  fzref2=repmat(fzref2,[1,size(fzref,2),1,1,1]);
       % fzref=fzref./fzref2;
        
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
        
    function val = readPar(fid_dir,field)
%val = readPar(fid_dir or parfile,field)



    
 if exist([fid_dir,'.fid'],'dir') 
    fid_dir = [fid_dir,'.fid'];
 elseif exist([fid_dir,'.img'],'dir') 
    fid_dir = [fid_dir,'.img'];
 end

 if ~isdir(fid_dir) && exist(fullfile(pwd,fid_dir),'file')
    fname=fid_dir;
 else
    
 if exist(fullfile(fid_dir,'procpar.orig'),'file')
  fname=fullfile(fid_dir,'procpar.orig');
 else
   fname=fullfile(fid_dir,'procpar');
 end
end

%a=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','delimiter',' ');
%a=textread(parfile,'%s%s%s%s%s%s%s%s%s%s%s',);
a=textread(fname,'%s');

ind = strmatch(field,a,'exact');
if length(ind) ~=1
    val=[];
    warning('Field not found');
    return;
end

type = a{ind+1};

ind = ind+11;
na = str2num(a{ind});
if isempty(na)
    error('Number of values error');
end

val = a{ind+1};

         if ~isempty(str2num(val))
            val = zeros(1,na);
              for i=1:na
               val(i) = str2num(a{ind+i});
              end
              
              if type == '6'
                  val = val/1000000;
              end
         elseif na==1
             val = a{ind+1};
         else
             val=a(ind+1:ind+na)';
         end
            
function k=fov_shift_kspace(k,dx,fov)
%k=fov_shift_kspace(k,dx,fov)
% x and fov are 2D vectors
sz=size(k);
if numel(sz)==2 
    sz(3)=1;
end

n=size(k,1);
p1=exp(1i*2*pi*dx(1)/fov(1)*(1:n)');
p1=repmat(p1,[1,sz(2:end)]);

n=size(k,2);
p2=exp(1i*2*pi*dx(2)/fov(2)*(1:n));
p2=repmat(p2,[sz(1),1,sz(3:end)]);

if length(dx)==3

n=size(k,3);
p3=exp(1i*2*pi*dx(3)/fov(3)*(1:n));
p3=shiftdim(p3,-1);
p3=repmat(p3,[sz(1:2),1,sz(4:end)]);

k=k.*p1.*p2.*p3;
else
    
k=k.*p1.*p2;
end

    
        
        