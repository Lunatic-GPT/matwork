function ManualAllineate_afni_SAVEALL(master_set,input,out_prefix,info_input)
% ignore the oblique information;
% save the number of intersecting voxels for each tissue type (i.e. roi
% value) in the input data.

% NN: take the maximum value of all values in input that inset with the
% voxel in master_set; the input values are assumed to be integers.

[d_master,info_master]=BrikLoad(master_set);

mat=zeros(3,4);
mat(:,1:3)=eye(3);

if ~exist('info_input','var')
[d,info_input]=BrikLoad(input);
else
    
    if isa(input,'char')
        d=ri(input);
    else
        d=input;
    end
end

n_out = uint8(zeros([size(d_master),max(d(:))]));

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

  

d=double(d);
inv_ijk2real_m=ijk2real_m(:,1:3)\eye(3);
inv_mat=mat(:,1:3)\eye(3);

mat3prod=inv_ijk2real_m*inv_mat*ijk2real;
mat3prod_b=-inv_ijk2real_m*inv_mat*mat(:,4);
mat2prod=-inv_ijk2real_m*ijk2real_m(:,4);

mat_3prod_2prod=mat3prod_b+mat2prod;

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
            
            if any(ijk_m<1) || any(ijk_m'>size(d_master))
                continue;
            end

            %d_out(ijk_m(1),ijk_m(2),ijk_m(3))=d_out(ijk_m(1),ijk_m(2),ijk_m(3))+d(i,j,k);
            
            if d(i,j,k)>0
              n_out(ijk_m(1),ijk_m(2),ijk_m(3),d(i,j,k))=n_out(ijk_m(1),ijk_m(2),ijk_m(3),d(i,j,k))+1;
            end
        end
    end
end

%d_out(n_out>0)=d_out(n_out>0)./n_out(n_out>0);


%{
opt.Prefix=prefix;

opt.OverWrite='y';
WriteBrik(d_out,info_master,opt);
%}

save(out_prefix,'n_out');


function pvec=find_permute_matrix(orientation)

pvec=zeros(1,3);
orientation=orientation';
 ind=find(orientation(:)=='L');
 pvec(ceil(ind/2))=1;
 
 ind=find(orientation(:)=='P');
 pvec(ceil(ind/2))=2;
    
 ind=find(orientation(:)=='I');
 pvec(ceil(ind/2))=3;
    
    
    
    
    
    
    


    
    
