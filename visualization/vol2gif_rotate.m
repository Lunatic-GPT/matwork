function vol2gif_rotate(d,dim,nstep,delayTime,prefix)
%% dim: 1*2; the first is the dimension for rotation; the second is the dimension for projection.

if ~exist('delayTime','var')
    delayTime = 0.2;
end

sz=size(d);
sz(dim(2))=[];
d2=zeros([sz,nstep]);


for i=1:nstep
    disp(i);
    tmp=rotate_images(d,dim(1),360/nstep*(i-1),true);    
    d2(:,:,i)=project(tmp,dim(2));
end

d3=uint8(d2/max(d2(:))*255);
vol2gif(d3,gray(255),prefix,delayTime);


function res=project(m,dim)

 alpha=0.9;
 m2=cumsum(m,dim);
    
    m3=alpha.^(m2-1).*m;
    
    res=squeeze(sum(m3,dim));
    
    
    