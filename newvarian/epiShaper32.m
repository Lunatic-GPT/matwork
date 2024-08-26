function epiShaper32(fid_dir,triple_ref,format,nophase)
% epiShaper32(fid_dir,triple_ref,format,nophase)
% nophase: do not save phase image.default true.
% format: 'a' or 's', default 'a';
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

orient=readPar(fid_dir,'orient');
orient=orient(2:end-1);

b=z2;
switch orient 
    case 'trans90'  %1,2,3 in b - l2r,u2d in stimulate - l2r, a2p,s2i in vnmrj
   
  % b=flipdim(b,2);  
   b=flipdim(b,1);  
  order = {'L2R','A2P','S2I'};
   
  case 'trans'
   b = permute(b,[2,1,3,4]);   %do I really need this? yes 
   b=flipdim(b,2); 
   b=flipdim(b,1);
   order={'L2R','A2P','S2I'};
   
 case 'sag90'  % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
  %  b = permute(b,[2,1,3,4]);
  
      b=flipdim(b,1);    
  %    b=flipdim(b,2);
      order = {'A2P','S2I','L2R'};
      
    case 'sag' % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
      b=permute(b,[2,1,3,4]);   
      b=flipdim(b,2); 
      b=flipdim(b,1);
      order = {'A2P','S2I','L2R'};
      
    case 'cor90'  
      b=  flipdim(b,1);
   %   b=flipdim(b,2);
      order = {'L2R','S2I','P2A'};
    case 'cor'
      b=permute(b,[2,1,3,4]);
     b=flipdim(b,2);
       b=flipdim(b,1);
      order = {'L2R','S2I','P2A'};
      
    otherwise
        
end

[sgn,pm]=reorient_sdt2afni(order);

tmp=b;
for i=1:3
    if sgn(i)==-1
        tmp=flipdim(tmp,i);
    end
end

afni=permute(tmp,[pm,4]);
    
if strcmp(format,'s')
 writesdt4(abs(b),[prefix,'_mag']);
 if ~nophase
  writesdt4(angle(b),[prefix,'_ph']);
 end
 return;
end

%%


thk = readPar(fid_dir,'thk');%in mm
tr = readPar(fid_dir,'tr');
pss=readPar(fid_dir,'pss'); %in cm
%img = slice_reorder(img,pss);    

    pro=readPar(fid_dir,'pro'); %in cm
    ppe=readPar(fid_dir,'ppe'); %in cm
    sz=size(afni);
    if length(pss)>1
        pss_sort = sort(pss);
        thk = (pss_sort(2)-pss_sort(1))*10;  %pss in cm 
    end
    
    
       
    %varian convention:
    %              axial;  axial90;  coronal; coronal90; sag;  sag90   
    %- to +    pro: A2P;    L2R;      S2I;     L2R;      S2I;  A2P;
    %- to +    ppe: R2L;    A2P;      R2L;     S2I;      P2A;  S2I;
    %- to +    pss: S2I;    S2I;      P2A;     P2A;      L2R;  L2R;
    %  
    %     
    % varian lab coordinates: [x,y,z] = [R2L,P2A,S2I] in roi panel;
    [sgn,pm]=reorient_varian2afni(orient);
    center = [pro,ppe,mean(pss)]*10;
    delta = [lro*10/nread,lpe*10/nphase,thk];
  
    center= center.*sgn;
    center=center(pm);
    
    delta=delta(pm);
    delta([1,3])=-delta([1,3]);
    if length(sz)==2
        sz(3)=1;
    end
    %    
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
    info.ORIGIN = center-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    info.DELTA = delta;
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I  % afni convention + is R2L, A2P, I2S
    write_afni(abs(afni),[prefix,'_mag'],info);
    if ~nophase
     write_afni(angle(afni),[prefix,'_ph'],info);
    end
    
    
function [sgn,pm]=reorient_varian2afni(orient)
        
        sgn=zeros(1,3);
        pm=zeros(1,3);
 
   
            
    switch orient
        case 'trans'
            order = {'A2P','R2L','S2I'};
        case 'trans90'
            order = {'L2R','A2P','S2I'};
        case 'cor'
            order = {'S2I','R2L','P2A'};
        case 'cor90'
            order = {'L2R','S2I','P2A'};
        case 'sag'
            order = {'S2I','P2A','L2R'};
        case 'sag90'
            order = {'A2P','S2I','L2R'};
        otherwise
            error('unknown orient');
    end
   
    %              axial;  axial90;  coronal; coronal90; sag;  sag90   
    %- to +    pro: A2P;    L2R;      S2I;     L2R;      S2I;  A2P;
    %- to +    ppe: R2L;    A2P;      R2L;     S2I;      P2A;  S2I;
    %- to +    pss: S2I;    S2I;      P2A;     P2A;      L2R;  L2R;
    
       for i=1:3 
        switch order{i}
            case 'R2L'
             sgn(i)=1;
             pm(1)=i; 
            case 'L2R'
              sgn(i)=-1;
              pm(1)=i;
            case 'A2P'
              sgn(i)=1;
              pm(2)=i;
            case 'P2A'
               sgn(i)=-1;
               pm(2)=i;
            case 'I2S'
               sgn(i)=1;
               pm(3)=i;
            case 'S2I'
              sgn(i)=-1;
              pm(3)=i;
            otherwise
                error('unknow directin');
                
        end
       end
       
       
       
     function [sgn,pm]=reorient_sdt2afni(orient)
        
        sgn=zeros(1,3);
        pm=zeros(1,3);
 
   if ~iscell(orient)
            
    switch orient
        case 'trans'
            order = {'A2P','R2L','S2I'};
        case 'trans90'
            order = {'L2R','A2P','S2I'};
        case 'cor'
            order = {'S2I','R2L','P2A'};
        case 'cor90'
            order = {'L2R','S2I','P2A'};
        case 'sag'
            order = {'S2I','P2A','L2R'};
        case 'sag90'
            order = {'A2P','S2I','L2R'};
        otherwise
            error('unknown orient');
    end
   
   else
       order=orient;
   end
    %L2R,A2P,S2I
    
       for i=1:3 
        switch order{i}
            case 'L2R'
             sgn(i)=1;
             pm(1)=i; 
            case 'R2L'
              sgn(i)=-1;
              pm(1)=i;
            case 'A2P'
              sgn(i)=1;
              pm(2)=i;
            case 'P2A'
               sgn(i)=-1;
               pm(2)=i;
            case 'S2I'
               sgn(i)=1;
               pm(3)=i;
            case 'I2S'
              sgn(i)=-1;
              pm(3)=i;
            otherwise
                error('unknow directin');
                
        end
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
        
        
        
        