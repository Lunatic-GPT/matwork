function mask=cstable_gems_mTR_LS(N,nTRcs,reduction,ns,nc,suffix)
%cstable_gems_mTR_LS(N,nTRcs,reduction,ns,nc,dorand,suffix)
% N: number of phase encoding lines (before reduction)
% nTRcs: number of TRs (paradigm TR factor*number of TR_per_trial)
% reduction: reduction factor
% ns: number of slices; always assume ns =1 for the moment.
% nc: number of fully sampled central k-space lines 
% suffix: a suffix to the phase table file name 
%{
N=64;
nTRcs=240;
reduction = 4;
ns=1;

nc=8;  %includes the edge line
%}
%%
dorand = true;
if mod(nTRcs,reduction)~=0
    nTRcs=reduction*ceil(nTRcs/reduction);
    warning('mod(nTRcs,reduction)~=0, set nTRcs to %d',nTRcs);
end

nout =(N/reduction-nc);
npe=zeros(nout+nc,ns,nTRcs);

oint = [1:N/2-nc/2,N/2+nc/2+1:N];

mask = zeros(N,ns,nTRcs);
  
for isl=1:ns
    
    % generate out k-space lines
    
     
for j=1:nTRcs/reduction
   k=randperm(length(oint));
     
    for i=1:reduction
      
        ind=N/2-nc/2+1:N/2+nc/2;
        tmp=k((i-1)*nout+1:i*nout);
        ind(end+1:end+nout)=oint(tmp);
        ind=sort(ind);
        npe(:,isl,(j-1)*reduction+i)=ind;
        %mask(:,ind)=1;
        mask(ind,isl,(j-1)*reduction+i)=1;
    
        
    end
       
end
%figure;imshow(squeeze(mask(:,isl,:)),[]);

end
%%
if dorand
npe2=zeros(size(npe));

for i=1:size(npe,2)
    for j=1:size(npe,3)
        tmp=randperm(size(npe,1));
      npe2(:,i,j)=npe(tmp,i,j);
    end
end
else
    npe2=npe;
end


fname = sprintf('petable_R%d_N%d_nTR%d_ns%d_nc%d_%s',reduction,N,nTRcs,ns,nc,suffix);

if dorand
    fname = [fname,'_rand'];
end

if exist(fname,'file')
    warning('%s exists',fname);
else
  fid=fopen(fname,'w');
  npe2=squeeze(npe2);
  npe2=npe2';
  npe2=npe2(:);
  %fprintf(fid,'t1 =\n');
  for i=1:length(npe2(:))
    fprintf(fid,'%d\n',npe2(i)-N/2-1);
  end
  fclose(fid);
  fprintf('petable saved in %s\n',fname);
end
