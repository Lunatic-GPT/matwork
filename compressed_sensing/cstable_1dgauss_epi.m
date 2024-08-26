function mask=cstable_1dgauss_epi(npe,reduction,ns,nTR,nc,fname_suffix)
% cstable_1dguass(npe,reduction,ns,nTR,nc[,fname_suffix])
% generate random k-space lines based on a gaussian distribution.
% always acquire k=0 line.
% npe: number of phase encoding lines before reduction
% reduction: reduction factor
% ns: number of slice
% nTR: number of TR
% fname_suffix: a string attached to the end of the default file name.
sig=0.5;
scl=fsolve(@(x) area_gauss(x,sig,npe)-1/reduction,0.5);
disp(scl);

petable2=zeros(npe/reduction,nTR,ns);
mask=zeros(npe,nTR,ns);

bound=floor((nc+1)/2);
for isl=1:ns
    
    i=0;
while i<nTR
    
petable=-npe/2;

while length(petable)<(npe/reduction-nc)
    
    if bound ==0
      x=(-npe/2:npe/2-1)/npe*2;
    else
      x=[-npe/2:-bound,bound:npe/2-1]/npe*2;
    end
        
    mtmp=rand(1,length(x))<scl.*exp(-x.^2/2/sig/sig);
  %  mask=rand(1,npe)<1/reduction;
    tmp=x(mtmp)/2*npe;
    tmp2=setdiff(tmp,petable);
    
    rind2=randperm(length(tmp2));
    petable=[petable,tmp2(rind2)];
    
end

petable=[-bound+1:bound-1,petable(1:npe/reduction-nc)];
petable=sort(petable);

if any(diff(petable)>10)
    continue;
end
i=i+1;
petable2(:,i,isl)=petable;

%disp(petable);
mask(petable+npe/2+1,i,isl)=1;
%mask=repmat(mask,[npe,1]);

end

figure;plot(sum(mask(:,:,isl),2));
figure;imshow(mask(:,:,isl),[]);
fprintf('Percent of lines unsampled during any TR = %f\n',100*(length(find(sum(mask(:,:,isl),2)==0))/npe));
end

if exist('fname_suffix','var') && ~isempty(fname_suffix) 
  fname=sprintf('blipfactor_1dgauss_N%d_R%d_ns%d_nTR%d_nc%d%s',npe,reduction,ns,nTR,nc,fname_suffix);
else
   fname=sprintf('blipfactor_1dgauss_N%d_R%d_ns%d_nTR%d_nc%d',npe,reduction,ns,nTR,nc); 
end

if exist(fname,'file')
    warning('%s already exist',fname);
else
 fid=fopen(fname,'w');
 blipfactor=petable2(2:end,:,:)-petable2(1:end-1,:,:);
 for i=1:length(blipfactor(:))
    fprintf(fid,'%d\n',blipfactor(i));
 end
 fclose(fid);
end

function frac=area_gauss(scl,sig,n)

a=(-n/2:n/2-1)/n*2;
scl=scl*ones(1,n);
scl(n/2+1)=1;
frac=sum(exp(-a.^2/sig/sig/2).*scl)/n;