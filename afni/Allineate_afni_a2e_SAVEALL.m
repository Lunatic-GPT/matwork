function Allineate_afni_a2e_SAVEALL(master_set,mat1D,input,prefix,info_input)
% save the number of intersecting voxels for each tissue type (i.e. roi
% value) in the input data.

% Use T1_al_e2a_only_mat.aff12.1D; for plumb orientation this file is not
% present; maybe same as the following file?
% In afni, use T1_al_mat.aff12.1D
% NN: take the maximum value of all values in input that inset with the
% voxel in master_set; the input values are assumed to be integers.


[d_master,info_master]=BrikLoad(master_set);

fid=fopen(mat1D,'r');
d=textscan(fid,'%f','CommentStyle','#');
mat=d{1};
fclose(fid);

mat=reshape(mat,[4,3])';

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

n_out = uint8(zeros([size(d_master),max(d(:))]));

  ijk2real=info_input.IJK_TO_DICOM_REAL;
  ijk2real=reshape(ijk2real,[4,3])';

ijk2real_m=info_master.IJK_TO_DICOM_REAL;
ijk2real_m=reshape(ijk2real_m,[4,3])';
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

save(prefix,'n_out');


function pvec=find_permute_matrix(orientation)

pvec=zeros(1,3);
orientation=orientation';
 ind=find(orientation(:)=='L');
 pvec(ceil(ind/2))=1;
 
 ind=find(orientation(:)=='P');
 pvec(ceil(ind/2))=2;
    
 ind=find(orientation(:)=='I');
 pvec(ceil(ind/2))=3;
    
    
    
    
    
    
    


    
    
