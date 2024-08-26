function  getDiffPar( fname,outf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mat=readbPar([fname,'/method'],'PVM_DwGradVec');
b=readbPar([fname,'/method'],'PVM_DwEffBval');
dir=mat./repmat(sqrt(sum(mat.^2,1)),[3,1]);
dir(isnan(dir))=0;
fid=fopen(outf,'w');
fprintf(fid,'b=\n');
for i=1:length(b)
    fprintf(fid, '%f   ',b(i));
end

fprintf(fid,'\n');

fprintf(fid,'dir=\n');

for j=1:3
for i=1:length(b)
    fprintf(fid,'%f   ',dir(j,i));
end
fprintf(fid,'\n');
end




fclose(fid);

