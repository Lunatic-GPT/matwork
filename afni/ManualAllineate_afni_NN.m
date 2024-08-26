function ManualAllineate_afni_NN(master_set,input,out_prefix,info_input)
% calculate the Ave of all voxels in input that intersect with the voxel in master_set
%ManualAllienate_afni_NN(master_set,input[,out_prefix,info_input])
if ~exist('out_prefix','var')
    out_prefix=[strtok(input,'+'),'_to_',strtok(master_set,'+')];
end

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

n_out = zeros([size(d_master),max(d(:))]);


xyz= Orientation2xyz(info_input.Orientation);
    
  ijk2real=zeros(3,4);
  for i=1:3
   ijk2real(i,xyz(i))=info_input.DELTA(xyz(i));
  end
  ijk2real(:,4)=info_input.ORIGIN(xyz);

  
xyz_m= Orientation2xyz(info_master.Orientation);
    
  ijk2real_m=zeros(3,4);
  for i=1:3
   ijk2real_m(i,xyz_m(i))=info_master.DELTA(xyz_m(i));
  end
  ijk2real_m(:,4)=info_master.ORIGIN(xyz_m);



%%
 


d=double(d);
inv_ijk2real_m=ijk2real_m(:,1:3)\eye(3);
inv_mat=mat(:,1:3)\eye(3);

mat3prod=inv_ijk2real_m*inv_mat*ijk2real;
mat3prod_b=-inv_ijk2real_m*inv_mat*mat(:,4);
mat2prod=-inv_ijk2real_m*ijk2real_m(:,4);

mat_3prod_2prod=mat3prod_b+mat2prod;

d_out=d_master(:,:,:,1)*0;

d=double(d);
for i=1:max(d(:))
    
    pos=ind2subb(size(d),find(d(:)==i))-1;
    
  for j=1:size(pos,1)
      
      sh=linspace(-0.4,0.4,3);
      for ix=1:3      
          for iy=1:3
              for iz=1:3
                  pos2=ijk2real*[pos(j,:)+[sh(ix),sh(iy),sh(iz)],1]';
                  
                  
                  pos3=mat(:,1:3)\(pos2-mat(:,4));
                  
                  
                  %pmat=find_permute_matrix(info_master.Orientation);
                  %pos3(pmat)=pos3;
                  
                  ijk_m=round(ijk2real_m(:,1:3)\(pos3-ijk2real_m(:,4)))+1;
                  
                  if any(ijk_m<1) || any(ijk_m'>size(d_out))
                      continue;
                  end
                  
                  d_out(ijk_m(1),ijk_m(2),ijk_m(3))=i;
              end
          end
      end
  end
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
    
    
    
    
    
    
    


    
    
