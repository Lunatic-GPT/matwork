function Allineate_Inverse_Mat(mat1D)

mat=load(mat1D);

mat=reshape(mat,[4,3])';

matp=inv(mat(:,1:3));
x0=mat(:,4);

x0p=-matp*x0;

mat_out=[matp,x0p]';

mat_out=mat_out(:)';

fid=fopen(['Inv_',mat1D],'w');

for i=1:length(mat_out)
    fprintf(fid,'%f   ',mat_out(i));
end
fclose(fid);

