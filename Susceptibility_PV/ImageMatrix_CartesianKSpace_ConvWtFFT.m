function res=ImageMatrix_CartesianKSpace_ConvWtFFT(mat,FOV,voxSize,voxSize_interp)
% res=Image_CartesianKSpace_ConvWtFFT(fr,center,FOV,voxSize,voxSize_interp)
% fr: the function of f(x,y) or f(x,y,z); 
% x,y,z are matrices with the same size. x, y, z are the first, second, and
% third dimensions, respectively.
% FOV: field of view; mm; the center of FOV is assumed to be [0,0,0];
% voxSize: voxSize for acquisition


nd=ndims(mat);

if length(FOV)==1
    FOV=FOV*ones(1,nd);
end

if length(voxSize)==1
    voxSize=voxSize*ones(1,nd);
end

if length(voxSize_interp)==1
    voxSize_interp=voxSize_interp*ones(1,nd);
end

nvox_mat=size(mat);

voxSize_mat=FOV./nvox_mat;  %voxSize of the non-smoothed data and output

nvox=round(nvox_mat.*voxSize_mat./voxSize/2)*2;
nvox_interp =round(nvox_mat.*voxSize_mat./voxSize_interp/2)*2;

nz=(nvox_mat-nvox)/2; % number of elements to set to zero on each side
nn=(nvox_mat-nvox_interp)/2; % number of elements to set to null on each side

fres0=fft3c(mat);

if nd==3
    fres0([1:nz(1),end-nz(1)+1:end],:,:)=0;
    fres0(:,[1:nz(2),end-nz(2)+1:end],:)=0;
    fres0(:,:,[1:nz(3),end-nz(3)+1:end])=0;
    
    fres0([1:nn(1),end-nn(1)+1:end],:,:)=[];
    fres0(:,[1:nn(2),end-nn(2)+1:end],:)=[];
    fres0(:,:,[1:nn(3),end-nn(3)+1:end])=[];
   
    res=ifft3c(fres0);
else
    fres0([1:nz(1),end-nz(1)+1:end],:)=0;
    fres0(:,[1:nz(2),end-nz(2)+1:end])=0;
    
    fres0([1:nn(1),end-nn(1)+1:end],:)=[];
    fres0(:,[1:nn(2),end-nn(2)+1:end])=[];
    
    res=ifft2c(fres0);
end




