function Allineate_afni_e2a_NN(master_set,mat1D,input,info_input)

% Use T1_al_e2a_only_mat.aff12.1D; for plumb orientation this file is not
% present; maybe same as the following file?
% In afni, use T1_al_mat.aff12.1D
% NN: take the maximum value of all values in input that inset with the
% voxel in master_set; the input values are assumed to be integers.
% input is the epi data;
% master_set is the T1 data



[d_master,info_master]=BrikLoad(master_set);

mat=load(mat1D);

mat=reshape(mat,[4,3])';

d_out=0*d_master;

if ~exist('info_input','var')
[d,info_input]=BrikLoad(input);
else
    d=ri(input);
end
  ijk2real=info_input.IJK_TO_DICOM_REAL;
  ijk2real=reshape(ijk2real,[4,3])';

%ijk2real(1,2:3)=0;
%ijk2real(2,1:2)=0;
%ijk2real(3,[1,3])=0;
info_master.BRICK_TYPES=3;
ijk2real_m=info_master.IJK_TO_DICOM_REAL;
ijk2real_m=reshape(ijk2real_m,[4,3])';
d=double(d);
mat_inv=mat(:,1:3)\eye(3);
ijk2real_inv=ijk2real(:,1:3)\eye(3);

for i=1:size(d_master,1)
    disp(i);
    for j=1:size(d_master,2)
        for k=1:size(d_master,3)
            
            
            pos=[i,j,k]-1;
            
            pos2=ijk2real_m*[pos,1]';
            
            
            pos3=mat_inv*(pos2-mat(:,4));
            
            %%{
            ijk=round(ijk2real_inv*(pos3-ijk2real(:,4)))+1;
            
            if any(ijk<1) || any(ijk'>size(d))
                continue;
            end
            
            d_out(i,j,k)=d(ijk(1),ijk(2),ijk(3));
%}
            
            
        end
        
    end
end

prefix=strtok(input,'+');
prefix=strtok(prefix,'.');

opt.Prefix=[prefix,'_ToT1'];

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
    
    
    
    
    
    
    


    
    
