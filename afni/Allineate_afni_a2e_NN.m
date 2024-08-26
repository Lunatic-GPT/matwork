function mask=Allineate_afni_a2e_NN(master_set,mat1D,input,info_input)

% Use T1_al_e2a_only_mat.aff12.1D; for plumb orientation this file is not
% present; maybe same as the following file?
% In afni, use T1_al_mat.aff12.1D
% NN: take the maximum value of all values in input that inset with the
% voxel in master_set; the input values are assumed to be integers.

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

ijk2real_m=info_master.IJK_TO_DICOM_REAL;
ijk2real_m=reshape(ijk2real_m,[4,3])';
d=double(d);
for i=1:max(d(:))
    
    pos=ind2subb(size(d),find(d(:)==i))-1;
    
  for j=1:size(pos,1)
    pos2=ijk2real*[pos(j,:),1]';

    %pmat=find_permute_matrix(info_input.Orientation);
    %pos2=pos2(pmat);
    
    
  % pos3=mat*[pos2',1]';
    
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

prefix=strtok(input,'+');
prefix=strtok(prefix,'.');

opt.Prefix=[prefix,'_myal'];

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
    
    
    
    
    
    
    


    
    
