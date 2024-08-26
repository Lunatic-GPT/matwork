function res = GRAPPA_new(kData,kCalib,kSize,lambda, iPAT)
% data and kCalib: pe1*pe2*coils*rep
% define sampling pattern manually


coils = size(kData,3); % get sizes

res = kData*0;
%[AtA] = corrMatrix(kCalib,kSize); % build coil correlation matrix
AtA = dat2AtA(kCalib, kSize); % build coil calibrating matrix (kSize(1)*kSize(2)*coils)*(kSize(1)*kSize(2)*coils)

for n=1:coils
  
    fprintf('.');
	tmp = ARC(kData, AtA,kSize, n,lambda,iPAT); % reconstruct single coil imag
    res(:,:,n,:)=tmp;
end
fprintf('\n');

function [res] = ARC(kData, AtA, kSize, c,lambda,iPAT)
[sx,sy,nCoil,rep] = size(kData);

n1=(sx+kSize(1)-1);

n2=(sy+kSize(2)-1);

nzp1=(kSize-1)/2;

kData = zpad(kData,[n1,n2,nCoil,rep]);

dummyK = zeros(kSize(1),kSize(2),nCoil);
dummyK((end+1)/2,(end+1)/2,c) = 1;
idxy = find(dummyK);  %the point in the pattern to estimate the value for


res = zeros(sx,sy,rep);


kData=permute(kData,[3,1,2,4]);   %change to coil*pe1*pe2*rep

%calculate the kernels

pat=getPatterns(kData(:,:,:,1),kSize,nCoil,iPAT); %kSize(1)*kSize(2)*coils*iPAT(1)*iPAT(2)

for y = 1:iPAT(2)
    for x=1:iPAT(1)
        
        if pat(nzp1(1)+1,nzp1(2)+1,1,x,y)>0           
            res(x:iPAT(1):end,y:iPAT(2):end,:)=kData(c,x+nzp1(1):iPAT(1):end-nzp1(1),y+nzp1(2):iPAT(2):end-nzp1(2),:);
            continue;
        end
    
        kernel = calibrate(AtA,kSize,nCoil,lambda,pat(:,:,:,x,y),idxy); %nPoints(1)*nPoints(2)

 
        ind=ind2subb(kSize,find(pat(:,:,1,x,y)));
        kData2=zeros(size(ind,1),nCoil,ceil(n1/iPAT(1)),ceil(n2/iPAT(2)),rep);
        
        for i=1:size(ind,1)
          xind=  x+ind(i,1)-1:iPAT(1):n1;
          yind=  y+ind(i,2)-1:iPAT(2):n2;
          kData2(i,:,1:length(xind),1:length(yind),:)=kData(:,xind,yind,:);    
        end
       
       
        sz=size(kData2);
        if length(sz)==4
            sz(5)=1;
        end
        kData3=reshape(kData2,[prod(sz(1:2)),prod(sz(3:5))]); %(nPoints(1)*nPoints(2)*nCoil)*(n1/iPAT(1)*n2/iPAT(2)*rep),
        tmp2=reshape(transpose(kernel(:))*kData3,sz(3:5));
        res(x:iPAT(1):end,y:iPAT(2):end,:)=tmp2(1:floor((sx-x)/iPAT(1))+1,1:floor((sy-y)/iPAT(2))+1,:);
   
    end
end


function pat=getPatterns(kData,kSize,nCoil,iPAT)

% pat is relative to the kData which is zero padded
nnz=0;

for y = 1:iPAT(2)
    for x=1:iPAT(1)  
        tmp=kData(1,x:iPAT(1):end,y:iPAT(2):end,1);
        if nnz<sum(abs(tmp(:))>0)
            nnz=sum(abs(tmp(:))>0);
            xs=x;
            ys=y;
        end        
    end
end

pat=zeros(kSize(1),kSize(2),nCoil,iPAT(1),iPAT(2));

for y = 1:iPAT(2)
    for x=1:iPAT(1)  

        indx=xs-x+1:iPAT(1):kSize(1);
        
        indy=ys-y+1:iPAT(2):kSize(2);
        
        indx(indx<=0)=[];
        indy(indy<=0)=[];
        
        pat(indx,indy,:,x,y)=1;
        
    end
end
        
        

function rawkernel = calibrate(AtA, kSize, nCoil,  lambda, sampling,idxy)
% AtA: (kSize(1)*kSize(2)*coils)*(kSize(1)*kSize(2)*coils)
% kernel: [kSize,nCoil]
if nargin < 6
	sampling = ones([kSize,nCoil]);
end


sampling(idxy) = 0;
idxA = find(sampling);

Aty = AtA(idxA,idxy); %(kSize(1)*kSize(2)*coils-1)*1
AtA = AtA(idxA,idxA); %(kSize(1)*kSize(2)*coils-1)*(kSize(1)*kSize(2)*coils-1)

%kernel = sampling*0;

lambda = norm(AtA,'fro')/size(AtA,1)*lambda;

rawkernel = inv(AtA + eye(size(AtA))*lambda)*Aty; %(kSize(1)*kSize(2)*coils-1)*1
%kernel(idxA) = rawkernel; 
%kernel = inv(AtA + eye(size(AtA))*lambda)*Aty; %(kSize(1)*kSize(2)*coils-1)*1

