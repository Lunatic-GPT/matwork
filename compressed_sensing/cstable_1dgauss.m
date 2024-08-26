function mask=cstable_1dgauss(npe,reduction,ns,nTR,sig,fname_suffix,centric)
% cstable_1dguass(npe,reduction,ns,nTR,sig[,fname_suffix])
% generate random k-space lines based on a gaussian distribution.
% always acquire k=0 line.
% npe: number of phase encoding lines before reduction
% reduction: reduction factor
% ns: number of slice
% nTR: number of TR
% nc: number of fully sampled center k-space lines
% fname_suffix: a string attached to the end of the default file name.
% remove nc; always set nc = 0; 6/15/2014.
% now can change sigma. x between [-1, 1];
%sig=0.5;


scl=fsolve(@(x) area_gauss(x,sig,npe)-1/reduction,0.5);
disp(scl);

nc=0;

if ~exist('centric','var')
    centric = false;
end

petable2=zeros(npe/reduction,nTR,ns);
mask=zeros(npe,nTR,ns);

bound=floor((nc+1)/2);
for isl=1:ns
for i=1:nTR
    
petable=[];
while length(petable)<(npe/reduction-nc)
    
    if bound ==0
      x=(-npe/2:npe/2-1)/npe*2;
    elseif mod(nc,2)==0
      x=[-npe/2:-bound-1,bound:npe/2-1]/npe*2;
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

if mod(nc,2)==0
petable=[-bound:bound-1,petable(1:npe/reduction-nc)];
else
   petable=[-bound+1:bound-1,petable(1:npe/reduction-nc)];
end

if ~centric
 petable=sort(petable);
else
 [xtmp,itmp]=sort(abs(petable),'descend');   
 petable=petable(itmp);
end

petable2(:,i,isl)=petable;
%disp(petable);
mask(petable+npe/2+1,i,isl)=1;
%mask=repmat(mask,[npe,1]);

end


%figure;plot(sum(mask(:,:,isl),2));
%figure;imshow(mask(:,:,isl),[]);
fprintf('Percent of lines unsampled during any TR = %f\n',100*(length(find(sum(mask(:,:,isl),2)==0))/npe));
end

if exist('fname_suffix','var') && ~isempty(fname_suffix) 
  fname=sprintf('cstable_1dgauss_N%d_R%d_ns%d_nTR%d_sig%3.2f%s',npe,reduction,ns,nTR,sig,fname_suffix);
else
   fname=sprintf('cstable_1dgauss_N%d_R%d_ns%d_nTR%d_sig%3.2f',npe,reduction,ns,nTR,sig); 
end

if exist(fname,'file')
    warning('%s already exist. Not overwritten',fname);
else
 fid=fopen(fname,'w');
 for i=1:length(petable2(:))
    fprintf(fid,'%d\n',petable2(i));
 end
 fclose(fid);
end

mask=permute(mask,[1,3,2]);
%save(fname,'mask');

function frac=area_gauss(scl,sig,n)

a=(-n/2:n/2-1)/n*2;
scl=scl*ones(1,n);
%scl(n/2+1)=1;
frac=sum(exp(-a.^2/sig/sig/2).*scl)/n;

