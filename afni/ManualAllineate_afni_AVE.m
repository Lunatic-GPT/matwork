function ManualAllineate_afni_AVE(master_set,input,out_prefix,info_input)
% calculate the Ave of all voxels in input that intersect with the voxel in master_set

[d_master,info_master]=BrikLoad(master_set);


mat=zeros(3,4);
mat(:,1:3)=eye(3);

%d_out=0*d_master;

if ~exist('info_input','var')
[d,info_input]=BrikLoad(input);
else
    
    if isa(input,'char')
        d=ri(input);
    else
        d=input;
    end
end

n_out = zeros(size(d_master));


xyz= Orientation2xyz(info_input.Orientation);
    
  ijk2real=zeros(3,4);
  for i=1:3
   ijk2real(xyz(i),i)=info_input.DELTA(i);
   ijk2real(xyz(i),4)=info_input.ORIGIN(i);
   
  end
  %ijk2real(:,4)=info_input.ORIGIN(xyz);

  
xyz_m= Orientation2xyz(info_master.Orientation);
    
  ijk2real_m=zeros(3,4);
  for i=1:3
   ijk2real_m(xyz_m(i),i)=info_master.DELTA(i);
   ijk2real_m(xyz_m(i),4)=info_master.ORIGIN(i);
   
  end
  %ijk2real_m(:,4)=info_master.ORIGIN(xyz_m);



%%
 


d=double(d);
inv_ijk2real_m=ijk2real_m(:,1:3)\eye(3);
inv_mat=mat(:,1:3)\eye(3);

mat3prod=inv_ijk2real_m*inv_mat*ijk2real;
mat3prod_b=-inv_ijk2real_m*inv_mat*mat(:,4);
mat2prod=-inv_ijk2real_m*ijk2real_m(:,4);

mat_3prod_2prod=mat3prod_b+mat2prod;

d_out=d_master(:,:,:,1)*0;
for i=1:size(d,1)
  
     disp(i);
    for j=1:size(d,2)
        for k=1:size(d,3)
    
            pos=[i,j,k]-1;

          %  pos2=ijk2real*[pos,1]';

          %  pos3=inv_mat*(pos2-mat(:,4));
          % ijk_m=round(inv_ijk2real_m*(pos3-ijk2real_m(:,4)))+1;

            
            ijk_m_raw=mat3prod*[pos,1]'+mat_3prod_2prod;
            ijk_m=round(ijk_m_raw)+1;
            
            if any(ijk_m<1) || any(ijk_m'>size(d_out))
                continue;
            end

            d_out(ijk_m(1),ijk_m(2),ijk_m(3))=d_out(ijk_m(1),ijk_m(2),ijk_m(3))+d(i,j,k);
            
     
              n_out(ijk_m(1),ijk_m(2),ijk_m(3))=n_out(ijk_m(1),ijk_m(2),ijk_m(3))+1;
     
        end
    end
end

d_out(n_out>0)=d_out(n_out>0)./n_out(n_out>0);


if ~exist('out_prefix','var')
    out_prefix=strtok(input,'+');
    out_prefix=[out_prefix,'rsTo',strtok(master_set,'+')];
end

opt.Prefix=out_prefix;

opt.OverWrite='y';
WriteBrik(d_out,info_master,opt);



function pvec=find_permute_matrix(orientation)

pvec=zeros(1,3);
orientation=orientation';
 ind=find(orientation(:)=='L');
 pvec(ceil(ind/2))=1;
 
 ind=find(orientation(:)=='P');
 pvec(ceil(ind/2))=2;
    
 ind=find(orientation(:)=='I');
 pvec(ceil(ind/2))=3;
    
    
    
    
    
    
    


    
    
