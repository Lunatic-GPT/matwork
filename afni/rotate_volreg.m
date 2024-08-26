function mask=rotate_volreg(master_info,mat1D,d,input_info)

% d should have the same data oreientation as output from BrikLoad

a=textscan(fid,'
mat=load(mat1D);

mat=reshape(mat,[4,3])';

d_out=0*d_master;

  ijk2real=input_info.IJK_TO_DICOM_REAL;
  ijk2real=reshape(ijk2real,[4,3])';


ijk2real_m=master_info.IJK_TO_DICOM_REAL;
ijk2real_m=reshape(ijk2real_m,[4,3])';
d=double(d);
for i=1:max(d(:))
    
    pos=ind2subb(size(d),find(d(:)==i))-1;
    
  for j=1:size(pos,1)
    pos2=ijk2real*[pos(j,:),1]';

    %pmat=find_permute_matrix(input_info.Orientation);
    %pos2=pos2(pmat);
    
    
  % pos3=mat*[pos2',1]';
    
    pos3=mat(:,1:3)\(pos2-mat(:,4));
    
    
    %pmat=find_permute_matrix(master_info.Orientation);
    %pos3(pmat)=pos3;
    
    ijk_m=round(ijk2real_m(:,1:3)\(pos3-ijk2real_m(:,4)))+1;
    
    if any(ijk_m<1) || any(ijk_m'>size(d_out))
        continue;
    end
    
    d_out(ijk_m(1),ijk_m(2),ijk_m(3))=i;
  end
end



function pvec=find_permute_matrix(orientation)

pvec=zeros(1,3);
orientation=orientation';
 ind=find(orientation(:)=='L');
 pvec(ceil(ind/2))=1;
 
 ind=find(orientation(:)=='P');
 pvec(ceil(ind/2))=2;
    
 ind=find(orientation(:)=='I');
 pvec(ceil(ind/2))=3;
    
    
    
    
    
    
    


    
    
