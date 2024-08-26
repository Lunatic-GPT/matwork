im=zeros(64,64);

im(20:44,28:37)=1;

im(20:23,28:30)=0;
figure;imshow(im,[]);
%% test shift
fim=fft2c(im);

ky=linspace(-0.5,0.5,65);
ky=repmat(ky(1:end-1),[64,1]);
kx=ky';
k2=[kx(:),ky(:),0*kx(:)];

k3=cat(4,kx,ky,0*kx);

shift=[0,1,0]*2*pi; % phase shift
shift=reshape(shift,[1,1,1,3]);

shift=repmat(shift,[size(fim),1,1]);
d3 = fim.*exp(-1i*sum(shift.*k3,4));  % 

im_shift=ifft2c(d3);

figure;imshow(abs(im_shift),[]);



%% test rotation
fim=fft2c(im);

ky=linspace(-0.5,0.5,65);
ky=repmat(ky(1:end-1),[64,1]);
kx=ky';
k2=[kx(:),ky(:),0*kx(:)];
    
    
%     ind=[3,1,2];
%     m = transform_matrix_rotation_arb_axis(90,0,10);
%     
    ind=[1,2,3];
    m = transform_matrix_rotation_arb_axis(0,0,10);   
%     
%     ind=[1,3,2];
%     m = transform_matrix_rotation_arb_axis(90,90,10); % rotating about y;
%     
    k2=k2(:,ind);
   k3=m*k2';
   k3=k3';
    
   matrix=[64,64,1];
nufft=NUFFT3D(k3,1,[0,0,0],matrix(ind),1,1);

im2=nufft'*fim(:);

figure;imshow(abs(squeeze(im2)),[]);





    
    