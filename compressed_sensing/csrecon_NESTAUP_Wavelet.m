function res = csrecon_NESTAUP_Wavelet(kData, lambda)

%kData is a 2D k-space data, unsampled points = 0.  Data will be normalized by max(abs(kData(:)));
    
kData=kData/max(abs(kData(:)));


opts = [];
opts.MaxIntIter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 100;


dataSize = size(kData);
imSize=dataSize;

op = Wavelet('Daubechies',4,6);
opts.U=@(x) WLtran(op,x,imSize);
opts.Ut=@(x) WLtranT(op,x,imSize);
opts.xplug=kData(:)*0;

m=abs(kData)>0;
A=@(x)  afun(x,m, imSize, 'normal');
At=@(x) afun(x,m,  imSize,'transp');

Ac = @(x) counter( A, x);
Atc= @(x) counter( At, x);
 
 opts.stoptest = 1; 
%opts.TypeMin='tv';
%opts.errFcn=@(x) errFcn(x,z,m,img0);
La=my_normest(A,At,length(kData(:)));

        muf=0.000001;
        res=NESTA_UP(Ac,Atc,kData(m>0),lambda,La,muf,opts);
   
    res=reshape(res,imSize);
    
    
   

function res=WLtran(op,x,imSize)

        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
%   imSize_dyd=[448,448];
x=reshape(x,imSize);
x2 = zpad(x,imSize_dyd);

res=op*x2;

res=res(:);

function res=WLtranT(op,x,imSize)

 imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
% imSize_dyd=[448,448];
x=reshape(x,imSize_dyd);

res=op'*x;

res=crop(res,imSize);

res=res(:);


function y = afun(x,m, imSize,  tflag)

    
   %     dataSize=size(m);
    if strcmp(tflag,'transp')
        y=0*m;
        y(m>0) = x;
        x = ifft2c(y);
       
        y = x(:);
    else
        
        x = reshape(x,imSize);

       % i1=ceil((imSize-dataSize)/2);
       % x=x(i1(1)+1:i1(1)+dataSize(1),i1(2)+1:i1(2)+dataSize(2));
        
        y = fft2c(x);
        y = y(m>0);
        
    end

