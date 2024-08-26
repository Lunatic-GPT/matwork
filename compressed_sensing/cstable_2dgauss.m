function mask=cstable_2dgauss(npe,reduction,nTR,nc,fname_suffix)
% cstable_2dguass(npe,reduction,nTR,nc[,fname_suffix])
% generate random k-space lines based on a gaussian distribution.
% npe: number of phase encoding lines before reduction:1*1 or 1*2
% reduction: reduction factor
% nTR: number of TR
% nc: number of center lines that is sampled.
% fname_suffix: a string attached to the end of the default file name.

if length(npe)==1
    npe(2)=npe;
end

if length(nc)==1
    nc(2)=nc;
end

sig=0.6;
scl=fsolve(@(x) area_gauss(x,sig,npe,nc)-1/reduction,0.5);
disp(scl);

petable2=zeros(npe(1)*npe(2)/reduction,nTR);
mask=zeros(npe(1)*npe(2),nTR);

dist=density_gauss(scl,sig,npe,nc);


for i=1:nTR
    
petable=[];

while length(petable)<(npe(1)*npe(2)/reduction-nc(1)*nc(2))
    tmp=find(rand(npe(1),npe(2))<dist&dist<1);
    tmp2=setdiff(tmp,petable);
    rind2=randperm(length(tmp2));
    petable=[petable;tmp2(rind2)];
end

petable=[find(dist==1),petable(1:(npe(1)*npe(2)/reduction-nc(1)*nc(2)))];
petable=sort(petable);
petable2(:,i)=petable;
%disp(petable);
mask(petable,i)=1;
%mask=repmat(mask,[npe,1]);

end

mask=reshape(mask,[npe(1),npe(2),nTR]);

figure;imshow(sum(mask(:,:,:)/nTR,3),[]);
figure;imshow(mask(:,:,1),[]);
fprintf('Percent of points unsampled during any TR = %f\n',100*(length(find(sum(mask,3)==0))/npe(1)/npe(2)));

if exist('fname_suffix','var') && ~isempty(fname_suffix) 
  fname=sprintf('cstable_2dgauss_N%d_%d_R%d_nTR%d_nc%d_%d%s',npe(1),npe(2),reduction,nTR,nc(1),nc(2),fname_suffix);
else
  fname=sprintf('cstable_2dgauss_N%d_%d_R%d_nTR%d_nc%d_%d',npe(1),npe(2),reduction,nTR,nc(1),nc(2));
end

if exist(fname,'file')
    warning('%s already exist.  File not overwritten',fname);
else
    
 fid=fopen(fname,'w');
 for i=1:length(petable2(:))
    fprintf(fid,'%d\n',petable2(i));
 end
 fclose(fid);
end

function frac=area_gauss(scl,sig,n,nc)

dist=density_gauss(scl,sig,n,nc);

frac=sum(dist(:))/length(dist(:));


function dist=density_gauss(scl,sig,n,nc)

bound=floor((nc+1)/2);

x1=(-n(1)/2:n(1)/2-1)/n(1)*2;

x2=(-n(2)/2:n(2)/2-1)/n(2)*2;

[x,y]=meshgrid(x2,x1);

dist=scl.*exp(-(x.^2+y.^2)/sig/sig/2);



dist(n(1)/2+1-bound(1)+1:n(1)/2+1+bound(1)-1,n(2)/2+1-bound(2)+1:n(2)/2+1+bound(2)-1)=1;









