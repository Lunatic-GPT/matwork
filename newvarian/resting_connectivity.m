function resting_connectivity(fid_prefix,tpattern,bmask,smask,tr,out_prefix)
% resting_connectivity(fid_prefix,tpattern,bmask,smask,out_prefix)
% smask is a mask for the seed region or 
% a 1*3 or a 2*3 matrix specifying the coordinates of the diagonal corners of a cube, 0 based;

%sdt2afni(fid_prefix);

%tr=readPar(fid_prefix,'tr');
%tr=2.5;

if ~exist('out_prefix','var')
    out_prefix='con';
end

if ~isempty(tpattern)
 cmd= sprintf('3dTshift -tpattern %s -rlt+ -TR %fs -prefix %s_tShft %s+orig',tpattern,tr,fid_prefix,fid_prefix);
 unix(cmd);
 fid_prefix=[fid_prefix,'_tShft'];
end

%cmd=sprintf('3dvolreg -prefix %s_tShft_vr %s+orig',fid_prefix,fid_prefix);
%unix(cmd);



if ~isempty(bmask) 
    
    
     
[a,info]=BrikLoad([fid_prefix,'+orig']);
 if ~strcmp(bmask,'all')
      
   m=BrikLoad(bmask);
   bm_label='bm';
 else
     sz=size(a);
    m=ones(sz(1:3));
    bm_label='all';
 end
out_prefix=[out_prefix,'_',bm_label];
g = mean_roi(a,m);
g=squeeze(g);
g = g-mean(g);
g=g(:);
save('resting_global.1D','g','-ascii');

cmd=sprintf('3dDetrend -vector resting_global.1D -prefix %s_%s_dtr -polort 3 %s+orig',fid_prefix,bm_label,fid_prefix);
%cmd=sprintf('3dDetrend -prefix %s_tShft_vr_dtr -polort 3 %s_tShft+orig',fid_prefix,fid_prefix);

unix(cmd);

b=BrikLoad(sprintf('%s_%s_dtr+orig',fid_prefix,bm_label));

else
    [b,info]=BrikLoad([fid_prefix,'+orig']);
end

fa=fft(b,[],4);

t=size(b,4)*tr;
np = ceil(0.2*t);
np2=floor(0.003*t);
fa(:,:,:,np+2:end-np)=0;
fa(:,:,:,2:np2)=0;
fa(:,:,:,end-np2+2:end)=0;

a2=ifft(fa,[],4);

if isa(smask,'char')
  sm=BrikLoad(smask);
else
    sz=size(a2);
   sm=zeros(sz(1:3));
  if size(smask,1)==1
     smask(2,:)=smask(1,:);
     smask=smask+1;
  end
  
  smask=sort(smask,1);
  sm(smask(1,1):smask(2,1),smask(1,2):smask(2,2),smask(1,3):smask(2,3))=1;
  
end

%disp(smask);

WriteBrikEZ(a2,info,'resting_connectivity',[fid_prefix,'_lowf'],'');
ref=mean_roi(a2,sm);
figure;plot(ref);
for i=1:3
if smask(1,i)==smask(2,i)
   out_prefix=sprintf('%s_%d',out_prefix,smask(1,i));
else
   out_prefix=sprintf('%s_%dt%d',out_prefix,smask(1,i),smask(2,i)); 
end
end


correlateAnalysis([fid_prefix,'_lowf+orig'],ref,out_prefix);


