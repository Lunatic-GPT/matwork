function fsems_recon(fid_prefix)
% tbgemsShaper(fid_dir,format)
% first dim: ro; second dim: pe; third dim: slice;
% 1/22/2011: Wrote it. Tested for single-slice 2D data. XZ
% format: output format: 'a': afni format. 's': sdt format.
np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
ns=readPar(fid_prefix,'ns');

pro=readPar(fid_prefix,'pro');
ppe=readPar(fid_prefix,'ppe');
lro=readPar(fid_prefix,'lro');
lpe=readPar(fid_prefix,'lpe');
%figure;imagesc(m2);drawnow;
arraydim=readPar(fid_prefix,'arraydim');
ne=arraydim;

z=read_fid([fid_prefix,'.fid']);

etl=readPar(fid_prefix,'etl');

z=reshape(z,[np/2,etl,ns,nv/etl,ne]);

z=permute(z,[1,2,4,3,5]);
z=reshape(z,[np/2,nv,ns,ne]);

pelist=readPar(fid_prefix,'pelist');

if ~isempty(pelist)
    z2=zeros(size(z));
    for i=1:length(pelist)
     z2(:,pelist(i)+64,:,:)=z(:,i,:,:);
    end
end

z=z2;
z=fov_shift_kspace(z,[0,-ppe],[lro,lpe]);
b=ifft2c(z);

d=fid_prefix;
orient = readPar(d,'orient');
orient=orient(2:end-1);

switch orient 
   case 'trans90'  %1,2,3 in b - l2r,u2d in stimulate - l2r, a2p,s2i in vnmrj
   
   b=flipdim(b,2);  
%   b=flipdim(b,1);  
  order = {'L2R','A2P','S2I'};
   
  case 'trans'
   b = permute(b,[2,1,3,4]);   %do I really need this? yes 
   b=flipdim(b,2); 
   order={'L2R','A2P','S2I'};
   
 case 'sag90'  % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
  %  b = permute(b,[2,1,3,4]);
      b=flipdim(b,1);    
      b=flipdim(b,2);
      order = {'A2P','S2I','L2R'};
      
    case 'sag' % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
      b=permute(b,[2,1,3,4]);   
      b=flipdim(b,2); 
      order = {'A2P','S2I','L2R'};
      
    case 'cor90'  
      b=  flipdim(b,1);
      b=flipdim(b,2);
      order = {'L2R','S2I','P2A'};
    case 'cor'
      b=permute(b,[2,1,3,4]);
      b=flipdim(b,2);
      
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

writesdt4(abs(b),[fid_prefix,'_mag']);
writesdt4(angle(b),[fid_prefix,'_ph']);


%%


lro = readPar(d,'lro');%in cm
lpe = readPar(d,'lpe');%in cm
thk = readPar(d,'thk');%in mm
tr = readPar(d,'tr');
pss=readPar(d,'pss'); %in cm
%img = slice_reorder(img,pss);    

    pro=readPar(d,'pro'); %in cm
    ppe=readPar(d,'ppe'); %in cm
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
    np=readPar(d,'np');
    nread=np/2;
    nphase=readPar(d,'nv');
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
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I  % afni convention -2+ is R2L, A2P, I2S
    write_afni(abs(afni),[fid_prefix,'_mag'],info);
    write_afni(angle(afni),[fid_prefix,'_ph'],info);
   
    
    function [sgn,pm]=reorient_varian2afni(orient)
    % coordinate conversion.  
    % covert it to afni convention of -2+ is R2L, A2P, I2S
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
    % data conversion
    % orient is the direction of the original data
    % pm is the matrix used to convert it to afni order of [1,3,5] or L2R,A2P,S2I
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
                