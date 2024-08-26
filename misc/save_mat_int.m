function save_mat_int(r,fname)
% save_mat_int(r,fname)
% save a 2 dimensional interger matrix in text file

fid = fopen(fname,'w');

for i=1:size(r,1)
    fprintf(fid,'%d ',r(i,:));
    fprintf(fid,'\n');
end
fclose(fid);