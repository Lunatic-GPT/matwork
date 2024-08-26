
sz=size(fd_170);
[k2,Data]=get_k_data_moco(reshape(fd_170,[sz(1:2),1,sz(3)]),[],[400,400,2],[1,2,3],[0,0,0]);
nufft3d=NUFFT3D(k2(:,1:3),1,[0,0,0],[sz(1:2),1],1,1);

im3d=nufft3d'*Data(:,:);

nufft=NUFFT(k2(:,1:2),1,[0,0],sz(1:2));

%res = NUFFT3D(k,w,shift,imSize,nblock,npool);

im=[];
for i=1:sz(3)
    tic;
im(:,:,i)=nufft'*Data(:,i);%rand(200,200);
toc;
end

figure;imshow(sos(im,3),[]);
