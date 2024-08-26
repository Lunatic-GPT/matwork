function a2=MotionArtifactSim(a,motion,voxsize,rot,ax)

%% motion[3,npe,npar] or [3,npar] (unit mm); in the second case, all npe are acquired in the same TR, such as SPACE or MPRAGE
% rot [npe,npar] or [npar]: rotation angles in degrees.  To reduce the
% calculation time, the rotation is assumed to occur only around one axis,
% given by ax [1,2,3];
% The first, second, and third dimensions are assumed to be ro, pe, and
% par.
% FOV in mm;

if ndims(motion)==2
   motion = reshape(motion,[3,1,size(motion,2)]);
   rot=repmat(rot(:)',[size(a,2),1]);
end

motion = motion./voxsize';  % in units of voxel size;
fa=fft2c(a);
fa=fft1c(fa,3);
fa=add_translation(fa,motion);

if any(abs(rot(:))>0)
    ktrue=kspace_coord(size(a),rot,ax);  % in units of radians
    
    a2=rotate_data(fa,ktrue,ax);
    
    a2=ifft1c(a2,ax);
    
else
    a2=ifft2c(fa);
    a2=ifft1c(a2,3);
    
end
%  function a2=resample_data(fa,ktrue,ax)
%         
%         a2=fa*0;
%     for i=1:size(fa,ax)
%         tic;
%         if ax==1
%             a2(i,:,:)= nufft2d(squeeze(ktrue(2:3,i,:,:)),squeeze(img(i,:,:)));
%         elseif ax==2
%             a2(:,i,:)= nufft2d(squeeze(ktrue([1,3],:,i,:)),squeeze(img(:,i,:)));
%         else
%             a2(:,:,i)= nufft2d(ktrue(1:2,:,:,i),img(:,:,i));
%         end
%         
%         fprintf('Time remaining = %d s\n',round(toc*(size(a,ax)-i)));
%         
%     end


 function a2=rotate_data(fa,ktrue,ax)
        
        a2=fa*0;
    for i=1:size(fa,ax)
        tic;
        if ax==1
            a2(i,:,:)= inufft2d(squeeze(ktrue(2:3,i,:,:)),squeeze(fa(i,:,:)));
        elseif ax==2
            a2(:,i,:)= inufft2d(squeeze(ktrue([1,3],:,i,:)),squeeze(fa(:,i,:)));
        else
            a2(:,:,i)= inufft2d(ktrue(1:2,:,:,i),fa(:,:,i));
        end
        
        fprintf('Time remaining = %d s\n',round(toc*(size(a,ax)-i)));
        
    end
    

function res=inufft2d(ktrue,data)
ktrue=reshape(ktrue,[2,size(ktrue,2)*size(ktrue,3)]);

st=nufft_init(ktrue',size(data),[5,5],round(size(data)*1.4));
res=nufft_adj(data(:),st);
res=res/sqrt(length(data(:)));

res=fftshift(res,1);
res=fftshift(res,2);


function res=nufft2d(ktrue,img)
ktrue=reshape(ktrue,[2,size(ktrue,2)*size(ktrue,3)]);

st=nufft_init(ktrue',size(img),[7,7],round(size(img)*2));
img=fftshift(img,1);
img=fftshift(img,2);

res=nufft(img,st)/sqrt(length(img(:)));
%res=res/sqrt(length(data(:)));
res=reshape(res,size(img));



function ktrue=kspace_coord(sz,rot,ax)
      
    % rot in degree;
kx=get_k_radian(sz(1));
ky=get_k_radian(sz(2));    
kz=get_k_radian(sz(3));    

rot=rot*pi/180;

ktrue=zeros([3,sz]);
ind=1:3;
ind(ax)=[];

for j=1:sz(2)
     for k=1:sz(3)   

          a=rot(j,k);
         
         rotmat=eye(3);
         rotmat(ind,ind)=[cos(a),sin(a);-sin(a),cos(a)];        
    
         ktrue(:,:,j,k)=rotmat*[kx(:),ky(j)*ones(length(kx),1),kz(k)*ones(length(kx),1)]';  
         
    end
end       

function res=add_translation(fa,motion)


kx=get_k_radian(size(fa,1));
ky=get_k_radian(size(fa,2));    
kz=get_k_radian(size(fa,3));
kz=reshape(kz,[1,1,length(kz)]);

phx=kx(:).*motion(1,:,:);
phy=ky.*motion(2,:,:);
phz=kz.*motion(3,:,:);

res=fa.*exp(1i*phx).*exp(1i*phy).*exp(1i*phz);





%% a is the image data (a 3D matrix)
%% motion is the stdev of motion (assuming gaussian distributio)(in units of voxel size) in the left-right and inferior-superior directions, assuming no motion in the anterior-posterior direction.
%% typical value 0.3.
%% This function assumes the first dimension of matrix "a" is along anterior-posterior.  Make sure this is the case!
%% output: the motion-corrupted image
%% rotation angle; ax - the axis of rotation
    function res=get_k_radian(l)
        
        res=-l/2:(l/2-1);
        
        res=res*2*pi/l;
        



