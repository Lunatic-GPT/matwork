function cropimage(imag,marg)
% [l,r,u,d]
data = imread(imag);
info = imfinfo(imag);


if exist('marg','var')
   
    l=marg(1);
    r=marg(2);
    u=marg(3);
    d=marg(4);
   data(:,end-r+1:end,:)=[]; 
   data(:,1:l,:)=[];
   data(end-d+1:end,:,:)=[]; 
   data(1:u,:,:)=[];
else
    
    thr=255;
    i=1;
    while ~any(data(:,i,:)<thr)    
        i=i+1;
    end
    data(:,1:i-1,:) = [];

    i=0;
    while ~any(data(:,end-i,:)<thr)
        i=i+1;
    end    
    data(:,end-i+1:end,:)=[];


    i=1;
    while ~any(data(i,:,:)<thr)    
        i=i+1;
    end
    data(1:i-1,:,:) = [];

    i=0;
    while ~any(data(end-i,:,:)<thr)
        i=i+1;
    end
    
    data(end-i+1:end,:,:)=[];
    
end
prefix=strtok(imag,'.');
imwrite(data,[prefix,'_crop.tif'],'Resolution',[info.XResolution,info.YResolution]);